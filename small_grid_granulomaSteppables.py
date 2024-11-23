from cc3d.core.PySteppables import *
import cc3d.CompuCellSetup as CompuCellSetup
import logging
import numpy as np
import random
import pandas as pd
import os



#logging.info("Simulation started")


heterogeneity_coeff = 1.0
#parameter changes are written in comments
class small_grid_granulomaSteppable(SteppableBasePy):

	def __init__(self, frequency=1):
		SteppableBasePy.__init__(self,frequency)
	def start(self):
		output_dir = self.output_dir
		if output_dir is not None:
			logfilepath = Path(output_dir).joinpath('simlog' + '.txt')
		logging.basicConfig(
			filename=logfilepath,  # Log file name
			level=logging.DEBUG,            # Set the logging level
			format='%(asctime)s - %(levelname)s - %(message)s',  # Log format
			filemode='w'  # 'w' to overwrite the file each time, 'a' to append
		)

 
# class mac_recruitment_Steppable(SteppableBasePy):
# 	def __init__(self, frequency=1):
# 		SteppableBasePy.__init__(self, frequency)
#
# 	def start(self):
# 		self.macro = self.cell_type.macro
# 		self.Tcell = self.cell_type.Tcell
#
# 	def step(self, mcs):
# 		try:
# 			chem_field = self.field.attr
# 			x = random.randint(90, 95)
# 			y = random.randint(13, 18)
# 			if 250<mcs and random.random()<0.1:
# 				cell = self.new_cell(self.macro)
# 				self.cell_field[x:x+3, y:y+3, 0] = cell
# 				if random.random()<heterogeneity_coeff:
# 					cell.targetSurface=15
# 					cell.lambdaSurface=10.0
# 				else:
# 					cell.targetSurface=20
# 					cell.lambdaSurface=10.0
#
# 			if mcs>1500 and mcs % 5==0 and random.random()<0.1:
# 				cell = self.new_cell(self.Tcell)
# 				self.cell_field[x:x+3, y:y+3, 0] = cell
# 		except Exception as e:
#             logging.error(f"small_grid_granulomaSteppableError occurred: {str(e)}")


class mac_recruitment_Steppable(SteppableBasePy):
	def __init__(self, frequency=1):
		SteppableBasePy.__init__(self, frequency)

	def start(self):
		# Initialize cell types
		self.macro = self.cell_type.macro
		self.Tcell = self.cell_type.Tcell
		# logging.info("mac_recruitment_Steppable started")

	def step(self, mcs):
		try:
			# Random recruitment of cells based on MCS (Monte Carlo Step)

			x = random.randint(80, 85)
			y = random.randint(13, 18)

			# logging.debug(f"Step {mcs}: Random recruitment at (x, y): ({x}, {y})")

			# Macro recruitment after MCS 250
			if 250 < mcs and random.random() < 0.1:
				cell = self.new_cell(self.macro)
				self.cell_field[x:x+3, y:y+3, 0] = cell
				# logging.info(f"Macro cell created at (x, y): ({x}, {y})")

				# Assign surface properties based on heterogeneity coefficient
				if random.random() < heterogeneity_coeff:
					cell.targetSurface = 15
					cell.lambdaSurface = 10.0
					# logging.debug(f"Cell surface set: targetSurface=15, lambdaSurface=10.0")
				else:
					cell.targetSurface = 20
					cell.lambdaSurface = 10.0
					# logging.debug(f"Cell surface set: targetSurface=20, lambdaSurface=10.0")

			# T-cell recruitment after MCS 1500
			if mcs > 1500 and mcs % 5 == 0 and random.random() < 0.1:
				cell = self.new_cell(self.Tcell)
				self.cell_field[x:x+3, y:y+3, 0] = cell
				# logging.info(f"T-cell created at (x, y): ({x}, {y})")
			logging.info(f"mac_recruitment_Steppable at {mcs}")
		except Exception as e:
			logging.error(f"Error in mac_recruitment_Steppable at MCS {mcs}: {str(e)}")
			
			
class SurfaceSteppable(SteppableBasePy):
	def __init__(self, frequency=1): #frequency from 1 to 200
		SteppableBasePy.__init__(self, frequency)
	def start(self):
		for cell in self.cell_list_by_type(self.MTB):	
			cell.targetSurface=5
			cell.lambdaSurface=50.0	
		for cell in self.cell_list_by_type(self.TCELL):	
			cell.targetSurface=12
			cell.lambdaSurface=10.0
		
	def step(self, mcs):
		if mcs>1500:
			for cell in self.cell_list_by_type(self.MTB):	
				cell.targetSurface=5
				cell.lambdaSurface=50.0	
			for cell in self.cell_list_by_type(self.TCELL):	
				cell.targetSurface=12
				cell.lambdaSurface=10.0


#Steppeble for killing infected macrophages by T cells. 	
df4=[]
class inf_mac_killingSteppable(SteppableBasePy):
	def __init__(self, frequency=500): #frequency from 1 to 200
		SteppableBasePy.__init__(self, frequency) 
	def step(self, mcs):
		try:
			type_map = {
					1: "mtb",
					2: "macro",
					3: "inf_mac",
					4: "Tcell",
					5: "m2"
				}
			cells_to_delete = []
			for cell in self.cell_list:
				if cell.type == self.TCELL:
					for neighbor, commonSurfaceArea in self.get_cell_neighbor_data_list(cell):
						if neighbor and random.random()<0.005:
							if neighbor.type == self.INF_MAC or neighbor.type == self.M2:
								kill_data = [cell.id, cell.xCOM, cell.yCOM, neighbor.id, type_map.get(neighbor.type), neighbor.xCOM, neighbor.yCOM, mcs]
								df4.append(kill_data)
								self.delete_cell(neighbor)
			logging.info(f"inf_mac_killingSteppable at {mcs}")
		except Exception as e:
			logging.error(f"inf_mac_killingSteppableError occurred: {str(e)}")
		dataframe4 = pd.DataFrame(df4, columns=['tcellid', 'tcellx', 'tcelly', 'neighborid','neighbor_type', 'neighborx', 'neighbory', 'mcs'])
		output_dir = self.output_dir
		if output_dir is not None:
			output_path4 = Path(output_dir).joinpath('cell_killing' + '.csv')
			with open(output_path4, 'w') as fout:
				dataframe4.to_csv(output_path4, index=False)


		
# To 	create two types of infected cells inf_mac and m2 with a certain probability. 						
class InfectionSteppable(SteppableBasePy):
	def __init__(self, frequency=1):
		SteppableBasePy.__init__(self, frequency)
	def step(self, mcs):
		try:
			cells_to_delete = []
			for cell in self.cell_list:
				if cell.type == self.MACRO:
					for neighbor, commonSurfaceArea in self.get_cell_neighbor_data_list(cell):
						if neighbor and neighbor.type == self.MTB:
							if random.random() < 0.75:
								cell.type = self.INF_MAC
								cell.dict['timer'] = 0
								cell.dict['timer_threshold'] = 1500
								if random.random()<heterogeneity_coeff:
									cell.targetSurface=15
									cell.lambdaSurface=10.0
								else:
									cell.targetSurface=20
									cell.lambdaSurface=10.0
								self.delete_cell(neighbor)
							else:
								cell.type = self.M2
								cell.dict['timer'] = 0
								cell.dict['timer_threshold'] = 1500
								if random.random()<heterogeneity_coeff:
									cell.targetSurface=15
									cell.lambdaSurface=10.0
								else:
									cell.targetSurface=20
									cell.lambdaSurface=10.0
								self.delete_cell(neighbor)
			logging.info(f"InfectionSteppable at {mcs}")
		except Exception as e:
			logging.error(f"InfectionSteppableError occurred: {str(e)}")
			
class CellTypeConverterSteppable(SteppableBasePy):
	def __init__(self, frequency=1):
		SteppableBasePy.__init__(self, frequency)
	def start(self):
		pass
	def get_cell_center(self, cell):
		return np.array([cell.xCOM, cell.yCOM, cell.zCOM])
	def step(self, mcs):
		try:
			cells_to_convert = []
			cells_to_delete = []
			for cell in self.cell_list_by_type(self.INF_MAC):
				cell.dict['timer'] = cell.dict.get('timer', 0) + 1
				if cell.dict['timer'] >= cell.dict['timer_threshold']:
					cells_to_convert.append(cell)
			for cell in cells_to_convert:
				random_iter = random.random()
				if random_iter < 0.1:
					cell.type = self.MACRO
					cell.dict['timer'] = 0
				else:
					center = self.get_cell_center(cell)
					self.delete_cell(cell)
					for bac_num in range(6):
						new_cell = self.new_cell(self.MTB)
						offset = np.random.uniform(-1, 1, 3)
						x = int(center[0] + offset[0])
						y = int(center[1] + offset[1])
						z = int(center[2] + offset[2])
						x = max(0, min(x, self.dim.x - 1))
						y = max(0, min(y, self.dim.y - 1))
						z = max(0, min(z, self.dim.z - 1))
						self.cell_field[x, y, z] = new_cell
			for cell in self.cell_list_by_type(self.M2):
				cell.dict['timer'] = cell.dict.get('timer', 0) + 1
				if cell.dict['timer'] >= cell.dict['timer_threshold']:
					cells_to_delete.append(cell)
			for cell in cells_to_delete:
				self.delete_cell(cell)
			logging.info(f"CellTypeConverterSteppable at {mcs}")
		except Exception as e:
			logging.error(f"CellTypeConverterSteppableError occurred: {str(e)}")

class Hypoxia_deathSteppable(SteppableBasePy):
	def __init__(self, frequency=1):
		SteppableBasePy.__init__(self, frequency)
	def get_cell_center(self, cell):
		return np.array([cell.xCOM, cell.yCOM, cell.zCOM])
	def step(self, mcs):
		try:
			oxy_limit = 25.0
			cell_to_convert = []
			chemical_field = self.field.oxy
			if mcs %10 == 0 :
				for cell in self.cell_list:
					if cell.type == self.MACRO or cell.type == self.TCELL or cell.type == self.M2:

						x, y, z = int(cell.xCOM), int(cell.yCOM), int(cell.zCOM)
						oxy_concentration = chemical_field[x, y, z]
						if oxy_concentration < oxy_limit and random.random() < 0.05:
							self.delete_cell(cell)
					elif cell.type == self.INF_MAC:

						x, y, z = int(cell.xCOM), int(cell.yCOM), int(cell.zCOM)
						oxy_concentration = chemical_field[x, y, z]
						if oxy_concentration < oxy_limit and random.random() < 0.05:
							cell_to_convert.append(cell)
				for cell in cell_to_convert:
					center = self.get_cell_center(cell)
					self.delete_cell(cell)
					for bac_num in range(2):
						new_cell = self.new_cell(self.MTB)
						offset = np.random.uniform(-1, 1, 3)
						x = int(center[0] + offset[0])
						y = int(center[1] + offset[1])
						z = int(center[2] + offset[2])
						x = max(0, min(x, self.dim.x - 1))
						y = max(0, min(y, self.dim.y - 1))
						z = max(0, min(z, self.dim.z - 1))
						self.cell_field[x, y, z] = new_cell
				logging.info(f"Hypoxia_deathSteppable at {mcs}")
		except Exception as e:
			logging.error(f"Hypoxia_deathSteppableError occurred: {str(e)}")
					
df1=[]
df2=[]
df3=[]
class OutputFileSteppable(SteppableBasePy):
	def __init__ (self, frequency=1):
		SteppableBasePy.__init__(self, frequency)
		
	def step(self, mcs):
		try:
			chemical_field1 = self.field.attr
			chemical_field2 = self.field.oxy
			type_map = {
					1: "mtb",
					2: "macro",
					3: "inf_mac",
					4: "Tcell",
					5: "m2"
				}
			# logging.info(f"OutputFileSteppable1 at {mcs}: {str(e)}")
			if mcs %50 ==0:
				for cell in self.cell_list:
					com = (cell.xCOM, cell.yCOM, cell.zCOM)
					attr_concentration = chemical_field1[com]
					oxy_concentration = chemical_field2[com]
					cell_data = [cell.id, type_map.get(cell.type, "unknown"), cell.xCOM, cell.yCOM, attr_concentration, oxy_concentration, mcs]
					df1.append(cell_data)
				dataframe1 = pd.DataFrame(df1, columns=['cell_id', 'cell_type', 'cell_xcoordinate','cell_ycoordinate', 'attr_concentration', 'oxy_concentration', 'mcs'])
			# logging.info(f"OutputFileSteppable2 at {mcs}: {str(e)}")
			if mcs %50 ==0:
				cols = 	np.arange(0.0,100.0, 5)
				oxy_con_list = []
				oxy_con_list.append(mcs)
				for i in np.arange(5.0,100.0, 5):
					oxy_con = chemical_field2[i, 15, 0]
					oxy_con_list.append(oxy_con)
				df2.append(oxy_con_list)
				dataframe2 = pd.DataFrame(df2, columns = cols)
			# logging.info(f"OutputFileSteppable3 at {mcs}: {str(e)}")
			if mcs %50 ==0:
				cols = 	np.arange(0.0,100.0, 5)
				attr_con_list = []
				attr_con_list.append(mcs)
				for i in np.arange(5.0,100.0, 5):
					attr_con = chemical_field1[i, 15, 0]
					attr_con_list.append(attr_con)
				df3.append(attr_con_list)
				dataframe3 = pd.DataFrame(df3, columns = cols)
			# logging.info(f"OutputFileSteppable4 at {mcs}: {str(e)}")
			output_dir = self.output_dir
			if output_dir is not None:

				output_path1 = Path(output_dir).joinpath('CellBasedOutput' + '.csv')
				with open(output_path1, 'w') as fout:
					dataframe1.to_csv(output_path1, index=False)

				output_path2 = Path(output_dir).joinpath('OxygenConcentration' + '.csv')
				with open(output_path2, 'w') as fout:
					dataframe2.to_csv(output_path2, index=False)

				output_path3 = Path(output_dir).joinpath('AttrConcentration' + '.csv')
				with open(output_path3, 'w') as fout:
					dataframe3.to_csv(output_path3, index=False)
			# logging.info(f"OutputFileSteppable5 at {mcs}")
		except Exception as e:
			logging.error(f"OutputFileSteppableError occurred: {str(e)}")
				
			
			
			
		
			
			
					
