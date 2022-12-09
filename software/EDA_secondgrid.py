import warnings
import os
import numpy as np
import pandas as pd
import scipy.optimize
import scipy.stats as st
from scipy.special import erf
import colorcet
import tqdm
import bokeh.io
from bokeh.plotting import figure, show
from bokeh.layouts import row
from bokeh.transform import jitter
import iqplot

data_path= "../datasets/"
df = pd.read_csv(os.path.join(data_path, 'caulobacter_growth_events.csv'), header=0)

# adding time since growth event beginning
gb=df.groupby(['bacterium', 'growth event'])
t =[gb.get_group(x) for x in gb.groups]

for i in t:
    start_time=i['time (min)'].values[0]
    i['growth_time']=i['time (min)']-start_time
df=pd.concat(t)
warnings.filterwarnings('once')
df

gb=df.groupby(['bacterium'])
first, second =[gb.get_group(x) for x in gb.groups]

plots = [
    bokeh.plotting.figure(
        x_axis_label="growth_time",
        y_axis_label="area (µm²)",
        frame_height=150,
        frame_width=200,
        title=f"growth event {gevent}",
    )
    for gevent in second["growth event"].unique()
]

for gevent, g in second.groupby("growth event"):
    plots[gevent].circle(
        x="growth_time",
        y="area (µm²)",
        source=g,
        size=2,
    )

bokeh.io.show(bokeh.layouts.gridplot(plots, ncols=5))


from IPython.core.display import HTML
HTML(myplot_html)






