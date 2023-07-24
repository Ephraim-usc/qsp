import math
import numpy as np
import pandas as pd

import unum
import unum.units as units


class System:
  def __init__(self, analytes):
    self.analytes = analytes
    self.compartments = list()
  
  def add_compartment(self, compartment, volumes):
    self.compartments.append(compartment)
  
  def add_flow(self, analyte, compartment_a, compartment_b, coefficient):
    self.flows.append()

  # rate: function that receives a dict of analyte concentrations, and returns the reaction rate
  # delta: dict of concentration changes of the analytes per reaction
  def add_reaction(self, compartment, rate, delta):
    self.reactions.append()
  
  def print(self):
    pass



def antigen_on(params):
  return 5 * params["adc"] * params["drug"]



