class Search():
  def __init__(self, names, values, limits = None, names_fixed = None, precision = 0.05):
    self.names = names
    self.values = dict(zip(names, values))
    if limits is None:
      limits = [(-math.inf, math.inf) for _ in names]
    self.limits = dict(zip(self.names, limits))
    if names_fixed is None:
      names_fixed = []
    self.names_fixed = names_fixed
    self.lrs = {name:0.1 for name in names}
    self.precision = precision
  
  def set(self, values):
    assert self.values.keys() == values.keys()
    self.values = values
  
  def limit(self, name):
    limit = self.limits[name]
    low = limit[0]; high = limit[1]
    if isinstance(low, str):
      low = eval(low, self.values.copy())
    if isinstance(high, str):
      high = eval(high, self.values.copy())
    return low, high
  
  def get(self):
    return self.values
  
  def get_up(self, name):
    value = self.values[name]
    lr = self.lrs[name]
    low, high = self.limit(name)
    if value < (low + high)/2:
      step = (value - low) * lr
    else:
      step = (high - value) * lr
    step = max(step, 1e-5)
    values = self.values.copy()
    values[name] = round(min(high - 1e-5, value + step), 5)
    return values
  
  def get_down(self, name):
    value = self.values[name]
    lr = self.lrs[name]
    low, high = self.limit(name)
    if value < (low + high)/2:
      step = (value - low) * lr
    else:
      step = (high - value) * lr
    step = max(step, 1e-5)
    values = self.values.copy()
    values[name] = round(max(low + 1e-5, value - step), 5)
    return values
  
  def faster(self, name):
    self.lrs[name] = min(0.5, self.lrs[name] * 1.5)
  
  def slower(self, name):
    self.lrs[name] = max(self.precision, self.lrs[name] * 0.5)
  
  def cold(self):
    for name in self.names:
      if (name not in self.names_fixed) and (self.lrs[name] > self.precision):
        return False
    return True



def maximize(func, search):
  x = search.get()
  y_max = func(**x)
  print(str(x) + " " + str(y_max), flush = True)
  
  x_prev = x.copy()
  for _ in range(100):
    for name in [name for name in search.names if name not in search.names_fixed]:
      x_up = search.get_up(name)
      y_up = func(**x_up)
      x_down = search.get_down(name)
      y_down = func(**x_down)
      
      if (y_up > max(y_down, y_max)):
        search.set(x_up)
        search.faster(name)
        y_max = y_up
      elif (y_down > max(y_up, y_max)):
        search.set(x_down)
        search.faster(name)
        y_max = y_down
      else:
        search.slower(name)
    
    x = search.get()
    print(str(x) + " " + str(y_max), flush = True)
    
    if x_prev == x and search.cold():
      break
    x_prev = x.copy()
  
  return search.get(), y_max
