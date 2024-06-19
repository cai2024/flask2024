import os
import glob
import importlib
import sys



package_name = 'demo.blueprint.modules'
files = glob.glob(os.path.dirname(os.path.abspath(__file__))  + '/pipeline*.py')

for file in files:    
    name = file[len(os.path.dirname(os.path.abspath(__file__)))+10:-3]
    module = importlib.import_module(".pipeline_"+name , package=package_name)
    globals()[name+"Form_default"] = getattr(module, name+"Form_default")
    globals()[name+"Form_custom"] = getattr(module, name+"Form_custom")


  

class PipelineForm:
    def __init__(self):
        script_dir = os.path.dirname(os.path.abspath(__file__))
        self.suffixes = [filename.split("pipeline_")[1][:-3] for filename in os.listdir(script_dir) if filename.startswith("pipeline_")]

    def create_instance(self, class_type):    
        cls = globals().get(class_type)
        return cls()

    





