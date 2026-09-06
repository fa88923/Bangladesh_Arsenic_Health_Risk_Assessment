import pandas as pd

df = pd.read_spss("hh.sav")
for col in df.columns:
    print(col)