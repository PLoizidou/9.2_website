#!/usr/bin/env python
# coding: utf-8

# # Homework 9.2: Caulobacter growth (65 pts)
# 
# [Data download](https://s3.amazonaws.com/bebi103.caltech.edu/data/caulobacter_growth_events.csv)
# 
# <hr />

# In this problem, we will study the growth and division of *Caulobacter crescentus* over time. The lab of [Norbert Scherer](http://schererlab.uchicago.edu) at the University of Chicago acquired these data and published the work in PNAS, which you can download [here](https://doi.org/10.1073/pnas.1403232111).
# 
# The clever experimental set-up allows imaging of single dividing cells in conditions that are identical through time. This is accomplished by taking advantage of a unique morphological feature of *Caulobacter*.  The mother cell is adherent to the a surface through its stalk. Upon division, one of the daughter cells does not have a stalk and is mobile. The system is part of a microfluidic device that gives a constant flow. So, every time a mother cell divides, the un-stalked daughter cell gets washed away. In such a way, the dividing cells are never in a crowded environment and the buffer is always fresh. This also allows for easier segmentation.
# 
# We define a growth event as a time period during which a bacterium grows before dividing. After a division, a new growth event begins. Based on image data kindly provided by Charlie Wright and [Sri Iyer-Biswas](https://iyerbiswas.com), I processed the images to get a set of growth events. The data are available here: [https://s3.amazonaws.com/bebi103.caltech.edu/data/caulobacter_growth_events.csv](https://s3.amazonaws.com/bebi103.caltech.edu/data/caulobacter_growth_events.csv).
# 
# It is well known that in ideal conditions, bacteria grow exponentially. One bacterium divides to form two, those two divide to form four, those four divide to form eight, and so on. But what about each individual growth event? How does a single cell grow?
# 
# 
# **a)** One theoretical model we will consider is that each the growth of an individual bacterium is linear. That is,
# 
# \begin{align}
# a(t) = a_0(1 + k t),
# \end{align}
# 
# where $a$ denotes the area observed in the microscope images.  An alternative model is that each individual bacterium grows exponentially, such that
# 
# \begin{align}
# a(t) = a_0\mathrm{e}^{kt}.
# \end{align}
# 
# Considering each growth event to be independent of all others, develop a generative model for each of the two theoretical models. Comment on any considerations you made and concerns you may have with your modeling procedures. 
# 
# **b)** Using this model, perform parameter estimates for $a_0$ and $k$ for each growth event separately. Think about how to display your results graphically and make informative graphics.
# 
# **c)** Compare the two models. Do you think growth is exponential or linear?

# ### Import libraries and prepare data

# In[228]:


import warnings
import os

import numpy as np
import pandas as pd
import scipy.optimize
import scipy.stats as st
from scipy.special import erf
import colorcet
import bebi103
import tqdm

import bokeh.io
from bokeh.plotting import figure, show
from bokeh.layouts import row
from bokeh.transform import jitter
import iqplot
bokeh.io.output_notebook()
get_ipython().run_line_magic('load_ext', 'blackcellmagic')


# In[17]:


data_path= "../Data/"
df = pd.read_csv(os.path.join(data_path, 'caulobacter_growth_events.csv'), header=0)
df


# In[18]:


# adding time since growth event beginning
gb=df.groupby(['bacterium', 'growth event'])
t =[gb.get_group(x) for x in gb.groups]

for i in t:
    start_time=i['time (min)'].values[0]
    i['growth_time']=i['time (min)']-start_time
df=pd.concat(t)
warnings.filterwarnings('ignore')
df


# In[67]:


gb=df.groupby(['bacterium'])
first, second =[gb.get_group(x) for x in gb.groups]
first


# In[68]:


gb1=first.groupby(['growth event'])
lst1 =[gb1.get_group(x) for x in gb1.groups]


# In[21]:


gb2=second.groupby(['growth event'])
lst2 =[gb2.get_group(x) for x in gb2.groups]


# ### EDA

# In[22]:


# Bacterium 1 Plot

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


# In[23]:


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


# In[24]:


f=iqplot.ecdf(first['area (µm²)'].values, title='Bacterium 1', conf_int=True)
s=iqplot.ecdf(second['area (µm²)'].values, title='Bacterium 2', conf_int=True)
bokeh.io.show(row(f,s))


# In[25]:


plots = [
    bokeh.plotting.figure(
        x_axis_label="growth_time",
        y_axis_label="area (µm²)",
        frame_height=150,
        frame_width=200,
        title=f"growth event {gevent}",
    )
    for gevent in first["growth event"].unique()
]

for gevent, g in first.groupby("growth event"):
    plots[gevent].circle(
        x="growth_time",
        y="area (µm²)",
        source=g,
        size=2,
    )

bokeh.io.show(bokeh.layouts.gridplot(plots, ncols=3))


# In[26]:


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

bokeh.io.show(bokeh.layouts.gridplot(plots, ncols=3))


# ## Part A : Considering each growth event to be independent of all others, develop a generative model for each of the two theoretical models. Comment on any considerations you made and concerns you may have with your modeling procedures.

# As many processes are involved in cell growth, we choose to model the area it occupies in the microscope as being Normally distributed. We also assume that the area of each cell we measure is independent of all of the other cells we measure, and further that the distribution we use to describe the cell area of any given egg is the same as any other. That is to say that the areas of the cells are independent and identically distributed, abbreviated i.i.d.

# MAth

# In[45]:


def gen_normal(mu, sigma, n):
    return  np.random.normal(mu, sigma, n)


# ### Functions for linear growth model

# In[48]:


def log_normal_linear(params, data):
    """
    data should contain 1) paramas and 2) both time and area. In this order.
    """
    t=data[0]
    area=data[1]

    a0, k, sigma = params

    mu = a0+a0*k*t

    if mu.any() < 0 or sigma <= 0 or a0<=0 or k<=0: # no parameter can be negative 
        return -np.inf
   
    return np.sum(st.norm.logpdf(area, mu, sigma))


def mle_iid_norm_linear(data):
    """Perform maximum likelihood estimates for parameters for i.i.d. norm measurements, parametrized by mu (ff, fb, koff, D) and sigma"""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        res = scipy.optimize.minimize(
            fun=lambda params, data: -log_normal_linear(params, data),
            x0=np.array([0.9, 0.8, 0.3]),   # 5 initial guesses in this order  ff, fb, koff, D, sigma
            args=(data,),
            method='Powell'
            # bounds=((0,1), (0,1), (0,np.inf),(0,np.inf), (0,np.inf))
        )
    if res.success:
        return res.x
    else:
        raise RuntimeError('Convergence failed with message', res.message)
        
def draw_parametric_bs_reps_mle_linear(
    mle_fun, gen_fun, data, args=(), reps=1, progress_bar=False
):
    """
    Draw parametric bootstrap replicates of maximum likelihood estimator.
    """

    t=data[0]
    area=data[1]
    a0, k, sigma = mle_fun(data)
    mu = a0+a0*k*t

    if progress_bar:
        iterator = tqdm.tqdm(range(reps))
    else:
        iterator = range(reps)
        
    return np.array([mle_fun([t,gen_fun(mu, sigma, n=len(mu))]) for _ in iterator])


# Bacterium 1 linear

# In[49]:


df_mle_normal_1_l = pd.DataFrame(
    columns=["bacterium", "model", "growth event","parameter","mle","lower","upper"]
)

params_bootstrap_1_l=[]

for i in lst1:
    data=[i['time (min)'][i['time (min)']>0], i['area (µm²)'][i['time (min)']>0]] 
    
    #obtaining the MLE of the  parameters 
    a0, k, sigma = mle_iid_norm_linear(data)
    
    #obtaining parametric bootstrap replicates for each experiment 
    fo_bs_reps_parametric_1_l = draw_parametric_bs_reps_mle_linear(
        mle_iid_norm_linear,
        gen_normal,
        data,
        args=(),
        reps=100,
        progress_bar=True,
    )
    params_bootstrap_1_l.append(fo_bs_reps_parametric_1_l)
    
    # obtaining the 95% CI upper and lower bound for each parameter
    [lower_a0, lower_k, lower_sigma], [
        upper_a0,
        upper_k,
        upper_sigma
    ] = np.percentile(fo_bs_reps_parametric_1_l, [2.5, 97.5], axis=0)
    
    sub_df_1_l = pd.DataFrame(
        {
            "bacterium": ["1", "1", "1"],
            "model": ["linear", "linear", "linear"],
            "growth event": str(i['growth event'].iloc[0]+1),
            "parameter": ["a0", "k", "σ"],
            "mle": [a0, k, sigma],
            "lower": [lower_a0, lower_k, lower_sigma],
            "upper": [upper_a0, upper_k, upper_sigma],
        }
    )
    df_mle_normal_1_l= pd.concat([df_mle_normal_1_l, sub_df_1_l])


# In[50]:


df_mle_normal_1_l


# Bacterium 2 linear 

# In[51]:


df_mle_normal_2_l = pd.DataFrame(
    columns=["bacterium", "model", "growth event","parameter","mle","lower","upper"]
)

params_bootstrap_2_l=[]

for i in lst2:
    data=[i['time (min)'][i['time (min)']>0], i['area (µm²)'][i['time (min)']>0]] 
    
    #obtaining the MLE of the  parameters 
    a0, k, sigma = mle_iid_norm_linear(data)
    
    #obtaining parametric bootstrap replicates for each experiment 
    fo_bs_reps_parametric_2_l = draw_parametric_bs_reps_mle_linear(
        mle_iid_norm_linear,
        gen_normal,
        data,
        args=(),
        reps=100,
        progress_bar=True,
    )
    params_bootstrap_2_l.append(fo_bs_reps_parametric_2_l)
    
    # obtaining the 95% CI upper and lower bound for each parameter
    [lower_a0, lower_k, lower_sigma], [
        upper_a0,
        upper_k,
        upper_sigma
    ] = np.percentile(fo_bs_reps_parametric_2_l, [2.5, 97.5], axis=0)
    
    sub_df_2_l = pd.DataFrame(
        {
            "bacterium": ["2", "2", "2"],
            "model": ["linear", "linear", "linear"],
            "growth event": str(i['growth event'].iloc[0]+1),
            "parameter": ["a0", "k", "σ"],
            "mle": [a0, k, sigma],
            "lower": [lower_a0, lower_k, lower_sigma],
            "upper": [upper_a0, upper_k, upper_sigma],
        }
    )
    df_mle_normal_2_l= pd.concat([df_mle_normal_2_l, sub_df_2_l])


# In[52]:


df_mle_normal_2_l


# In[53]:


df_mle_linear = pd.concat([df_mle_normal_1_l, df_mle_normal_2_l])
df_mle_linear


# ### Functions for Exponential growth model

# In[218]:


dataZ=[lst1[1]['time (min)'][lst1[1]['time (min)']>=0], lst1[1]['area (µm²)'][lst1[1]['time (min)']>=1]] #using only data from t=0 and onwards


# In[219]:


t=dataZ[0]
area=dataZ[1]


# In[220]:


def theor_area(a0, k, t):
    """Compute mu"""
    return a0*np.exp(k*t)


def log_normal_exp(params, t, area):
    """
    data should contain 1) paramas and 2) both time and area. In this order.
    """
    a0, k, sigma = params

    mu = theor_area(a0, k, t)

    if mu.any() < 0 or sigma <= 0 or a0<=0 or k<=0: # no parameter can be negative 
        return -np.inf
   
    return np.sum(st.norm.logpdf(area, mu, sigma))


def resid(params, t, area):
    """Residual."""
    return area - theor_area(*params, t)


def mle_iid_norm_lstsq_exp(dataZ):
    """Compute MLE for parameters."""
    # Unpack data
    t=dataZ[0]
    area=dataZ[1]

    res = scipy.optimize.least_squares(
        resid, np.array([780, 0.9]), args=(t, area)
    )

    # Compute residual sum of squares from optimal params
    rss_mle = np.sum(resid(res.x, t, area)**2)

    # Compute MLE for sigma
    sigma_mle = np.sqrt(rss_mle / len(t))

    return tuple([x for x in res.x] + [sigma_mle])


def draw_parametric_bs_reps_mle_exp(
    mle_fun, gen_fun, data, args=(), reps=1, progress_bar=False
):
    """
    Draw parametric bootstrap replicates of maximum likelihood estimator.
    """

    t=data[0]
    area=data[1]
    a0, k, sigma = mle_fun(data)
    mu = a0*np.exp(k*t)

    if progress_bar:
        iterator = tqdm.tqdm(range(reps))
    else:
        iterator = range(reps)
        
    return np.array([mle_fun([t,gen_fun(mu, sigma, n=len(mu))]) for _ in iterator])


# In[221]:


mle_iid_norm_lstsq_exp(dataZ)


# In[222]:


mle_iid_norm_lstsq_exp(dataZ)


# In[128]:


mle_iid_norm_lstsq_exp(dataZ)


# In[132]:


mle_iid_norm_lstsq_exp(dataZ)


# In[136]:


mle_iid_norm_lstsq_exp(dataZ)


# In[140]:


mle_iid_norm_lstsq_exp(dataZ)


# Bacterium 1 exp

# In[110]:


df_mle_normal_1_e = pd.DataFrame(
    columns=["bacterium", "model", "growth event","parameter","mle","lower","upper"]
)

params_bootstrap_1_e=[]

for i in lst2:
    data=[i['time (min)'][i['time (min)']>0], i['area (µm²)'][i['time (min)']>0]] 
    
    #obtaining the MLE of the  parameters 
    a0, k, sigma = mle_iid_norm_lstsq_exp(data)
    
    #obtaining parametric bootstrap replicates for each experiment 
    fo_bs_reps_parametric_1_e = draw_parametric_bs_reps_mle_exp(
        mle_iid_norm_lstsq_exp,
        gen_normal,
        data,
        args=(),
        reps=100,
        progress_bar=True,
    )
    params_bootstrap_1_e.append(fo_bs_reps_parametric_1_e)
    
    # obtaining the 95% CI upper and lower bound for each parameter
    [lower_a0, lower_k, lower_sigma], [
        upper_a0,
        upper_k,
        upper_sigma
    ] = np.percentile(fo_bs_reps_parametric_1_e, [2.5, 97.5], axis=0)
    
    sub_df_1_e = pd.DataFrame(
        {
            "bacterium": ["1", "1", "1"],
            "model": ["exponential", "exponential", "exponential"],
            "growth event": str(i['growth event'].iloc[0]+1),
            "parameter": ["a0", "k", "σ"],
            "mle": [a0, k, sigma],
            "lower": [lower_a0, lower_k, lower_sigma],
            "upper": [upper_a0, upper_k, upper_sigma],
        }
    )
    df_mle_normal_1_e= pd.concat([df_mle_normal_1_e, sub_df_1_e])


# In[226]:


a0, k, sigma = mle_iid_norm_linear(data)
t=data[0]
area=data[1]
mu_l = a0+a0*k*t
#a0, k, sigma = mle_iid_norm_exp(data)
#mu_e = a0*np.exp(k*t)
p = figure(title="Bacterium", x_axis_label='Time (min)', y_axis_label='Area')
p.line(t, mu_l, line_color=(1,1,1), line_width=5, legend_label='linear' )
#p.line(t, mu_e, line_color=(0,0,1), line_width=5, legend_label='exponential' )
p.line(t,area, legend_label='data', line_width=4, line_color=(255,20,147))
# p.legend.location = "bottom_right"
# plots.append(p)
# show(row(plots))
show(p)


# assessing linear model

# In[231]:


rg = np.random.default_rng()


def sample_frap(params, sigma, t, size=1):
    """Generate samples of mean intensity vs time."""
    samples = np.empty((size, len(t)))

    for i in range(size):
        mu = theor_area(a0, k, t)
        samples[i] = rg.normal(mu, sigma)

    return samples


# Theoretical time points for predictive curve
t_theor = np.concatenate((np.linspace(-5, -0.01, 50), np.linspace(0, 25, 200)))

# Build predictive regression plots one-by-one
plots = []
for trial in np.sort(first["growth event"].unique()):
    params = df_mle_linear.loc[trial, :].values
    sigma = params[-1]
    params = params[:-1]

    plots.append(
        bebi103.viz.predictive_regression(
            samples=sample_frap(params, sigma, t_theor, size=1000),
            samples_x=t_theor,
            data=second.loc[
                df["growth event"] == trial, ["time (min)", "area (µm²)"]
            ].values,
            x_axis_label="time (min)",
            y_axis_label="area (µm²)",
            frame_height=150,
            frame_width=200,
            x_range=[-5, 25],
            title=f"growth event {trial}",
        )
    )

bokeh.io.show(bokeh.layouts.gridplot(plots, ncols=3))


# ## Part C Compare the two models. Do you think growth is exponential or linear?

# In[102]:


def Akaike(l_e, l_l, data):
    """
    input: the 3 MLEs for each experiment and the data for each experiment ([[time][area]])
    returns Akaike weight for the exponential and linear model. in this order.
    """
    l_e = log_normal_exp(l_e, data)
    l_l = log_normal_linear(l_l, data)
    AIC1 = -2 * l_e + 6
    AIC2 = -2 * l_l + 6
    w_e = 1 / (1 + np.exp(-0.5 * (AIC2 - AIC1)))
    w_l = 1 / (1 + np.exp(-0.5 * (AIC1 - AIC2)))
    return w_e, w_l


# In[103]:


l_e=mle_iid_norm_exp(data)
l_l=mle_iid_norm_linear(data)
Akaike(l_e, l_l,data)


# The Akaike weights (for all 107 independent growth events of the two bacteria) strongly suggest that the exponential model is superior to the linear one.

# In[ ]:




