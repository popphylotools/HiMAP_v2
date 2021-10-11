#!/usr/bin/env python
import os
import configparser
config = configparser.ConfigParser()
from context import project_dir
config.read(os.path.join(project_dir, 'config.ini'))

with open("./himap/welcome.txt") as f: 
    print (f.read())

for path in [config['Paths'].get(k) for k in config['Paths'].keys()]:
    print(path)
    os.makedirs(os.path.join(project_dir, path), exist_ok=True)
