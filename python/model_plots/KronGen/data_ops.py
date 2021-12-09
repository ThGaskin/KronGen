import numpy as np
import pandas as pd
import csv

from utopya.plotting import register_operation

def save_data(data, out_path: str):

    dimension = ''
    for key in data.coords:
        if (data.coords[key].size > 1):
            dimension = key
    df = pd.DataFrame(data, index = data.coords[dimension].data)
    df.to_csv(out_path+"/"+data.name+"_data.csv")

register_operation(name="save_data", func=save_data)
