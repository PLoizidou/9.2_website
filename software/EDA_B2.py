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
warnings.filterwarnings('once')
df

gb=df.groupby(['bacterium'])
first, second =[gb.get_group(x) for x in gb.groups]

# Bacterium 2 Plot

p2 = bokeh.plotting.figure(
    x_axis_label="growth_time",
    y_axis_label="area (µm²)",
    frame_height=400,
    frame_width=450,
)

# Colors for trails
colors = colorcet.b_glasbey_category10

for gevent, g in second.groupby("growth event"):
    p2.circle(
        x="growth_time",
        y="area (µm²)",
        source=g,
        size=2,
        color=colors[gevent],
        legend_label=f"growth event {gevent}",
    )

p2.title.text = 'Bacterium 2'
p2.add_layout(p2.legend[0], 'right')
p2.legend.spacing = 1
p2.legend.padding = 1
p2.legend.click_policy="mute"

bokeh.io.show(p2)


# myplot_html = file_html(p, CDN)
# # this HTML code is very long (~30 K), the cell below doesn't show all the code in NBviewer
# print(myplot_html) 


# In[5]:


from IPython.core.display import HTML
HTML(myplot_html)






