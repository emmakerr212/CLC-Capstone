# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 10:16:00 2026

@author: d4mcbrid
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# Reduction stage sensitivity analysis
# looking at: outlet CO2 profiles over time and energy analysis from reduction
ts = np.load('Pilot_red_base\\times.npy')


Cs_base = np.load('Pilot_red_base\\Cs.npy')
Cs_flow_5L = np.load('Pilot_red_ Flow 0.5\\Cs.npy')
Cs_flow_75 = np.load('Pilot_red_ Flow 0.75\\Cs.npy')
Cs_flow_5H = np.load('Pilot_red_ Flow 1.5\\Cs.npy')
Cs_flow_25 = np.load('Pilot_red_ Flow 1.25\\Cs.npy')
youtsdry = np.load('Pilot_red_base\\youtsdry.npy')

plt.figure()
plt.plot(ts,youtsdry[0],'tab:red',label='CH$_4$')
plt.plot(ts,youtsdry[1],'tab:purple',label='H$_2$')
plt.plot(ts,youtsdry[2],'tab:green',label='CO')
plt.plot(ts,youtsdry[3],'tab:blue',label='CO$_2$')
plt.legend()
plt.ylabel('Outlet Fraction (dry basis)')
plt.xlabel('Time (s)')

plt.figure()
plt.plot(ts,Cs_flow_5L[3,-1,:],'tab:green',linestyle='dotted',label='0.029')#' kg/m$^2$/s')
plt.plot(ts,Cs_flow_75[3,-1,:],'tab:orange',linestyle='dashed',label='0.044')# kg/m$^2$/s')
plt.plot(ts,Cs_base[3,-1,:],'tab:blue',linestyle='solid',label='0.058')# kg/m$^2$/s')
plt.plot(ts,Cs_flow_25[3,-1,:],'tab:red',linestyle='dashdot',label='0.073')# kg/m$^2$/s')
plt.plot(ts,Cs_flow_5H[3,-1,:],'tab:purple',linestyle=(0, (3, 5, 1, 5, 1, 5)),label='0.087')# kg/m$^2$/s')
plt.xlabel('Time (s)')
plt.ylabel('CO$_2$ Concentration (mol/m$^3$)')
plt.legend(title='Inlet Fuel Flux (kg/m$^2$/s)',handlelength=2.5)

Cs_OC_3 = np.load('Pilot_red_ OC 0.3\\Cs.npy')
Cs_OC_4 = np.load('Pilot_red_ OC 0.4\\Cs.npy')
Cs_OC_6 = np.load('Pilot_red_ OC 0.6\\Cs.npy')
Cs_OC_7 = np.load('Pilot_red_ OC 0.7\\Cs.npy')

plt.figure()
plt.plot(ts,Cs_OC_3[3,-1,:],'tab:green',linestyle='dotted',label='0.3')# kg Fe$_2$O$_3$/kg Red Mud')
plt.plot(ts,Cs_OC_4[3,-1,:],'tab:orange',linestyle='dashed',label='0.4')# kg Fe$_2$O$_3$/kg Red Mud')
plt.plot(ts,Cs_base[3,-1,:],'tab:blue',linestyle='solid',label='0.5')# kg Fe$_2$O$_3$/kg Red Mud')
plt.plot(ts,Cs_OC_6[3,-1,:],'tab:red',linestyle='dashdot',label='0.6')# kg Fe$_2$O$_3$/kg Red Mud')
plt.plot(ts,Cs_OC_7[3,-1,:],'tab:purple',linestyle=(0, (3, 5, 1, 5, 1, 5)),label='0.7')# kg Fe$_2$O$_3$/kg Red Mud')
plt.xlabel('Time (s)')
plt.ylabel('CO$_2$ Concentration (mol/m$^3$)')
plt.legend(title='  OC Concentration\n(kg Fe$_2$O$_3$/kg Red Mud)',title_fontsize=11,handlelength=2.5)

Cs_Hfuel15 =  np.load('Pilot_red_ Fuel H\\Cs.npy')
Cs_Lfuel85 = np.load('Pilot_red_ Fuel L\\Cs.npy')
Cs_HHfuel = np.load('Pilot_red_ Fuel HH\\Cs.npy')
Cs_LLfuel = np.load('Pilot_red_ Fuel LL\\Cs.npy')


header = " "*5 + f"{'CH₄':^10}{'H₂':^10}{'CO':^10}{'    CO₂':^11}{'        N₂':^11}"

labels = [f"{0.06:^10.3f}{0.10:^10.2f}{0.23:^10.3f}{0.14:^11}{0.47:^11.3f}",
          f"{0.07:^10.3f}{0.12:^10.2f}{0.28:^10.3f}{0.14:^11}{0.39:^11.3f}",
          f"{0.085:^10.3f}{0.14:^10.2f}{0.335:^10.3f}{0.14:^11}{0.3:^11.3f}",
          f"{0.1:^10.3f}{0.16:^10.2f}{0.385:^10.3f}{0.14:^11}{0.215:^11.3f}",
          f"{0.11:^10.3f}{0.18:^10.2f}{0.44:^10.3f}{0.14:^11}{0.13:^11.3f}"]
handles = [Line2D([0],[0],color='tab:green', lw=2),
           Line2D([0],[0],color='tab:orange', lw=2),
           Line2D([0],[0],color='tab:blue', lw=2),
           Line2D([0],[0],color='tab:red', lw=2),
           Line2D([0],[0],color='tab:purple', lw=2),]


plt.figure()
plt.plot(ts,Cs_LLfuel[3,-1,:],'tab:green',linestyle='dotted',label='Case 1')
plt.plot(ts,Cs_Lfuel85[3,-1,:],'tab:orange',linestyle='dashed',label='Case 2')
plt.plot(ts,Cs_base[3,-1,:],'tab:blue',linestyle='solid',label='Baseline')
plt.plot(ts,Cs_Hfuel15[3,-1,:],'tab:red',linestyle='dashdot',label='Case 3')
plt.plot(ts,Cs_HHfuel[3,-1,:],'tab:purple',linestyle=(0, (3, 5, 1, 5, 1, 5)),label='Case 4')
plt.xlabel('Time (s)')
plt.ylabel('CO$_2$ Concentration (mol/m$^3$)')
plt.legend(handlelength=2.5)
#legend = plt.legend(handles, labels, title=header, loc='lower right',frameon=True)
#plt.text(800,3.25,'Fuel Composition',ha='center',va='center',fontsize=10)
#legend.get_title().set_fontfamily("monospace")

# creating the tornado plot

import pandas as pd
from highlight_text import ax_text

labels = np.char.array([
    'Inlet Mass Flux',#n -50% - +50%',
    'OC Concentration',#\n 0.3 - 0.7 kg Fe$_2$O$_3$/kg Red Mud',
    'Fuel Concentration'])#'\n Low - High'])
# midpoint = 56.014 #MJ

# low_values = np.array([23.731,57.866,37.473]) #MJ

# high_values = np.array([88.533,51.463,80.365]) #MJ

midpoint = 56.014/3.6 #kWh

low_values = np.array([23.731/3.6,57.866/3.6,37.473/3.6]) #kWh

high_values = np.array([88.533/3.6,51.463/3.6,80.365/3.6]) #kWh

var_effect = np.abs(high_values - low_values)/midpoint

data = pd.DataFrame({'Labels':labels,
                     'Low values':low_values,
                    'High values':high_values,
                    'Variable effect':var_effect
                    })
data = data.sort_values(
    'Variable effect',
    ascending=True,
    inplace=False,
    ignore_index=False,
    key=None)

def tornado_chart(labels, midpoint, low_values, high_values, title='<Low> vs <High> sensitivity'):
    color_low = 'tab:blue'
    color_high = 'tab:red'
    ys = range(len(data['Labels']))[::1]
    plt.figure()
    
    for y, low_value, high_value in zip(ys, low_values, high_values):
        low_width = midpoint - low_value
        high_width = high_value - midpoint
        
        plt.broken_barh([(low_value,low_width),(midpoint,high_width)], (y-0.35,0.7), facecolors=[color_low,color_high],
                        edgecolors=['black','black'], linewidth=0.5)
        offset = 3
        
        if high_value > low_value:
            x_high = midpoint + high_width + offset
            x_low = midpoint - low_width - offset
        else:
            x_high = midpoint + high_width - offset
            x_low = midpoint - low_width + offset
       # plt.text(x_high, y, str(round(high_value)), va='center', ha='center')
      #  plt.text(x_low, y, str(round(low_value)), va='center', ha = 'center')
    # for y, sub in zip(ys,['0.3-0.7 kg Fe$_2$O$_3$/kg Red Mud',' ','0.029 - 0.087 kg/m$^2$/s']):
    #     plt.text( min(low_values)-20,y-0.22,sub,fontsize=10,ha='center',va='center',color='gray')
    
    plt.axvline(midpoint,color='black',linewidth=1)
    
    plt.text(1,1,'Fuel Composition', ha='center',va='center',color='black',fontsize=12)
    plt.text(1,0,'OC Concentration', ha='center', va='center',color='black',fontsize=12)
    plt.text(1,2,'Inlet Mass Flux',ha='center',va='center',color='black',fontsize=12)
    #plt.text(1,2-0.22,'<>')
  #  plt.text(1,1-0.2,'CH$_4$: 0.06, H$_2$: 0.10, CO: 0.23, CO$_2$: 0.14, N$_2$: 0.47', ha='center', va='center', color='tab:blue',alpha=0.75)
  #  plt.text(1,1-0.4,'CH$_4$: 0.11, H$_2$: 0.18, CO: 0.44, CO$_2$: 0.14, N$_2$: 0.13', ha='center', va='center', color='tab:red',alpha=0.75)
    ax_text(1,1-0.2,'<-30%> - <+30%>', ha='center', va='center',color='black',fontsize=10,
            highlight_textprops=[{'color':'tab:blue','alpha':0.75},{'color':'tab:red','alpha':0.75}])
    ax_text(1,2-0.2, '<0.029> - <0.087> kg/m²/s',ha='center',va='center',color='black',fontsize=10,
            highlight_textprops=[{'color':'tab:blue','alpha':0.75},{'color':'tab:red','alpha':0.75}])
    ax_text(1,-0.2, '<0.3> - <0.7> kg Fe₂O₃/kg Red Mud',ha='center',va='center',color='black',fontsize=10,
            highlight_textprops=[{'color':'tab:blue','alpha':0.75},{'color':'tab:red','alpha':0.75}])
    ax = plt.gca()
    ax.spines[['right','left','top']].set_visible(False)
    ax.set_yticks([])
    
    ax_text(x = midpoint, y = len(labels), s = title, color='black',fontsize=15,va='center',ha='center',
            highlight_textprops=[{'color':color_low,'fontweight':'bold'},{'color':color_high,'fontweight':'bold'}],ax=ax)
    ax.spines['left'].set_position(('data', min(low_values) - 7))
    plt.text(56.014/3.6,2.65,'Baseline = 15.56',fontsize=10,ha='center',va='center')
    plt.xlabel('Energy (kWh)')
   # plt.yticks(ys,labels,fontsize=14,ha='right')
    plt.xlim(20*1.3/3.6,75*1.3/3.6)
    plt.ylim(-0.5,len(labels)-0.5)
    plt.tight_layout()
    plt.tick_params(left=False)
    plt.show()
    return

tornado_chart(data['Labels'],midpoint,data['Low values'],data['High values'])
