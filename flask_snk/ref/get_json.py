import json
import sys

# 假设命令行第一个参数是JSON文件名，第二个参数是输出文件名
json_file = sys.argv[1]
json_key = sys.argv[2]


key_list=json_key.split(',')

with open(json_file, 'r') as file:
    data = json.load(file)

    if len(key_list)==2:
        my_data = data[key_list[0]][key_list[1]]

    elif len(key_list)==3:
        my_data = data[key_list[0]][key_list[1]][key_list[2]]


print(my_data)
    
