import os
from .pipeline_WGS import WGSForm_default,WGSForm_custom
from .pipeline_CS import CSForm_default,CSForm_custom
#from .pipeline_WGBS import WGBSForm_default,WGBSForm_custom



class PipelineForm:
    class_map = {
        'WGS_default': WGSForm_default,
        'WGS_custom': WGSForm_custom,
        'CS_default': CSForm_default,
        'CS_custom': CSForm_custom,

        #'WGBS_default': WGBSForm_default,
        #'WGBS_custom': WGBSForm_custom,
       
    }
    def __init__(self):
        script_dir = os.path.dirname(os.path.abspath(__file__))
        self.suffixes = find_pipeline_suffixes(script_dir)

    def create_instance(self, class_type):   
            return self.class_map[class_type]()





def find_pipeline_suffixes(directory):
    suffixes = []
    for filename in os.listdir(directory):
        if filename.startswith("pipeline_"):
            suffix = filename.split("pipeline_")[1][:-3]  # 分割字符串获取后缀
            suffixes.append(suffix)
    return suffixes

