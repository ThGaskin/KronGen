import csv
import numpy as np
import pandas as pd

from utopya.plotting import register_operation

# Save data to a csv file
def save_data(data, out_path: str):

    # Find sweep dimension
    dimension = ''
    for key in data.coords:
        if (data.coords[key].size > 1):
            dimension = key

    # Create data frame and save to csv
    df = pd.DataFrame(data, index = data.coords[dimension].data)
    df.to_csv(out_path+"/"+data.name+"_data.csv")

register_operation(name="save_data", func=save_data)
