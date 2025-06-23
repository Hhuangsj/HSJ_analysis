from simtk import openmm, unit
from simtk.openmm import app
import numpy as np
class DualTemperatureLangevinIntegrator(openmm.CustomIntegrator):
    def __init__(self, temperature_cold, temperature_hot, hot_atoms, friction_coeff, step_size):
        """
        基于LangevinMiddleIntegrator的双温区积分器
        
        参数:
        temperature_cold: 冷区温度 (K)
        temperature_hot: 热区温度 (K)
        hot_atoms: 高温原子索引集合
        friction_coeff: 摩擦系数 (1/ps)
        step_size: 时间步长 (ps)
        """
        super().__init__(step_size)
        
        # 单位转换
        kB = unit.MOLAR_GAS_CONSTANT_R  # kJ/(mol·K)
        kT_cold = kB * temperature_cold
        kT_hot = kB * temperature_hot
        kT_cold_val = kT_cold.value_in_unit(unit.kilojoule_per_mole)
        kT_hot_val = kT_hot.value_in_unit(unit.kilojoule_per_mole)
        
        # 添加全局变量
        self.addGlobalVariable("a", 0)  # exp(-friction*dt)
        self.addGlobalVariable("b", 0)  # sqrt(1-exp(-2*friction*dt))
        self.addGlobalVariable("friction", friction_coeff)
        self.addGlobalVariable("dt", step_size)
        self.addGlobalVariable("kT_cold", kT_cold_val)
        self.addGlobalVariable("kT_hot", kT_hot_val)
        
        # 添加每原子温度变量
        self.addPerDofVariable("local_kT", kT_cold_val)
        self.addPerDofVariable("x1", 0)  # 位置缓存
        
        # 添加积分步骤
        self.addUpdateContextState()
        
        # v = v + dt*f/m
        self.addComputePerDof("v", "v + dt*f/m")
        self.addConstrainVelocities()
        
        # x = x + 0.5*dt*v
        self.addComputePerDof("x", "x + 0.5*dt*v")
        
        # 计算Langevin参数
        self.addComputeGlobal("a", "exp(-friction*dt)")
        self.addComputeGlobal("b", "sqrt(1-exp(-2*friction*dt))")
        
        # 应用Langevin热浴: v = a*v + b*sqrt(local_kT/m)*gaussian
        self.addComputePerDof("v", "a*v + b*sqrt(local_kT/m)*gaussian")
        
        # x = x + 0.5*dt*v
        self.addComputePerDof("x", "x + 0.5*dt*v")
        
        # 保存位置
        self.addComputePerDof("x1", "x")
        
        # 约束位置
        self.addConstrainPositions()
        
        # 修正速度: v = v + (x-x1)/dt
        self.addComputePerDof("v", "v + (x-x1)/dt")
        self.addConstrainVelocities()
        
        # 保存热区原子列表和热区温度值
        self.hot_atoms = hot_atoms
        self.kT_hot_val = kT_hot_val
        self.kT_cold_val = kT_cold_val

    def setHighTemperatureAtoms(self, context):
        """设置高温区域原子的局部温度"""
        # 获取local_kT变量的值 (3D向量)
        local_kT = self.getPerDofVariableByName("local_kT")
        
        # 创建新的local_kT数组
        for i in range(len(local_kT)):
            if i in self.hot_atoms:
                # 高温原子设置为热区温度值 (保留3个分量)
                local_kT[i] = (self.kT_hot_val, self.kT_hot_val, self.kT_hot_val)
            else:
                # 低温原子设置为冷区温度值 (保留3个分量)
                local_kT[i] = (self.kT_cold_val, self.kT_cold_val, self.kT_cold_val)
        
        # 设置新的local_kT值
        self.setPerDofVariableByName("local_kT", local_kT)

class DualTemperatureReporter:
    """自定义报告器，记录系统总温度及分区域温度"""
    def __init__(self, file, reportInterval, hot_atoms, total_atoms):
        """
        file: 输出文件
        reportInterval: 报告间隔步数
        hot_atoms: 高温原子索引列表
        total_atoms: 系统总原子数
        """
        self._reportInterval = reportInterval
        self._out = open(file, 'w') if isinstance(file, str) else file
        self._hot_atoms = set(hot_atoms)
        self._needsPositions = True
        self._needsVelocities = True
        self._needsForces = False
        self._needsEnergy = True
        
        # 设置原子信息
        self._total_atoms = total_atoms
        self._cold_atoms = set(range(total_atoms)) - self._hot_atoms
        self._n_cold = len(self._cold_atoms)
        self._n_hot = len(self._hot_atoms)
        
        # 写入表头
        headers = [
            "Step", "Time(ps)", 
            "Potential Energy(kJ/mol)", 
            "Kinetic Energy(kJ/mol)", 
            "Total Energy(kJ/mol)", 
            "Temperature(K)", 
            "Hot Region Temp(K)", 
            "Cold Region Temp(K)",
            "Speed(ns/day)"
        ]
        self._out.write("#" + ",".join(headers) + "\n")
    
    def describeNextReport(self, simulation):
        steps = self._reportInterval - simulation.currentStep % self._reportInterval
        return (steps, self._needsPositions, self._needsVelocities, 
                self._needsForces, self._needsEnergy)
    
    def report(self, simulation, state):
        """收集并报告模拟状态"""
        try:
            # 获取当前步骤信息
            step = simulation.currentStep
            time = state.getTime()
            
            # 获取能量信息
            energy = state.getKineticEnergy() + state.getPotentialEnergy()
            ke = state.getKineticEnergy()
            pe = state.getPotentialEnergy()
            
            # 获取速度信息 (as numpy数组)
            velocities = state.getVelocities(asNumpy=True)
            
            # 计算整体温度
            if self._total_atoms > 0:
                total_temp = (2*ke/(3*unit.MOLAR_GAS_CONSTANT_R*self._total_atoms)).value_in_unit(unit.kelvin)
            else:
                total_temp = 0.0
                
            # 计算分区温度
            kB = unit.MOLAR_GAS_CONSTANT_R  # kJ/(mol·K)
            hot_ke = 0.0 * unit.kilojoule_per_mole
            cold_ke = 0.0 * unit.kilojoule_per_mole
            
            # 计算高温区动能 (使用numpy索引)
            if self._n_hot > 0:
                hot_indices = list(self._hot_atoms)
                for i in hot_indices:
                    v = velocities[i]
                    mass = simulation.system.getParticleMass(i)
                    # 使用索引访问速度分量而不是属性
                    hot_ke += 0.5 * mass * (v[0]*v[0] + v[1]*v[1] + v[2]*v[2])
            
            # 计算低温区动能 (使用numpy索引)
            if self._n_cold > 0:
                cold_indices = list(self._cold_atoms)
                for i in cold_indices:
                    v = velocities[i]
                    mass = simulation.system.getParticleMass(i)
                    # 使用索引访问速度分量而不是属性
                    cold_ke += 0.5 * mass * (v[0]*v[0] + v[1]*v[1] + v[2]*v[2])
            
            # 计算分区温度 (公式: T = 2*KE/(3N*kB))
            if self._n_hot > 0:
                hot_temp = (2*hot_ke/(3*kB*self._n_hot)).value_in_unit(unit.kelvin)
            else:
                hot_temp = 0.0
                
            if self._n_cold > 0:
                cold_temp = (2*cold_ke/(3*kB*self._n_cold)).value_in_unit(unit.kelvin)
            else:
                cold_temp = 0.0
            
            # 计算运行速度
            if step > 0:
                speed = (time.value_in_unit(unit.picosecond) / step) * 86400 * 1e-6  # ns/day
            else:
                speed = 0
            
            # 格式化输出
            data = [
                step,
                time.value_in_unit(unit.picosecond),
                pe.value_in_unit(unit.kilojoule_per_mole),
                ke.value_in_unit(unit.kilojoule_per_mole),
                energy.value_in_unit(unit.kilojoule_per_mole),
                total_temp,
                hot_temp,
                cold_temp,
                speed
            ]
            
            self._out.write(",".join(map(str, data)) + "\n")
            self._out.flush()
            
        except Exception as e:
            print(f"Error in reporter: {str(e)}")
    
    def __del__(self):
        if hasattr(self, '_out') and self._out:
            self._out.close()
