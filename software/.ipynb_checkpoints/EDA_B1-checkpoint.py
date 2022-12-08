import warnings
import os
import numpy as np
import pandas as pd
import colorcet
import bokeh.io
from bokeh.plotting import figure, show
from bokeh.transform import jitter
from IPython.core.display import HTML
from bokeh.resources import CDN
from bokeh.embed import file_html

data_path= "../datasets/"
df = pd.read_csv(os.path.join(data_path, 'caulobacter_growth_events.csv'), header=0)

# adding time since growth event beginning
gb=df.groupby(['bacterium', 'growth event'])
t =[gb.get_group(x) for x in gb.groups]

for i in t:
    start_time=i['time (min)'].values[0]
    i['growth_time']=i['time (min)']-start_time
df=pd.concat(t)
warnings.filterwarnings('ignore')
df

gb=df.groupby(['bacterium'])
first, second =[gb.get_group(x) for x in gb.groups]

p = bokeh.plotting.figure(
    x_axis_label="growth_time",
    y_axis_label="area (µm²)",
    frame_height=400,
    frame_width=450,
)

# Colors for trails
colors = colorcet.b_glasbey_category10

for gevent, g in first.groupby("growth event"):
    p.circle(
        x="growth_time",
        y="area (µm²)",
        source=g,
        size=2,
        color=colors[gevent],
        legend_label=f"growth event {gevent}",
    )

p.title.text = 'Bacterium 1'
p.add_layout(p.legend[0], 'right')
p.legend.spacing = 1
p.legend.padding = 1
p.legend.click_policy="mute"

bokeh.io.show(p)




