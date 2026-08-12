# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 10:16:00 2026

@author: d4mcbrid
"""

import numpy as np
import matplotlib.pyplot as plt

# Oxidation stage sensitivity analysis
# looking at: temperature profiles and energy analysis from oxidation
ts = np.load('Pilot_ox_base\\times.npy')
zs = np.load('Pilot_ox_base\\length.npy')
zs = zs[3:-1]-zs[3]
Xs_base = np.load("Pilot_ox_base\\Xs.npy")


Ts_base = np.load('Pilot_ox_base\\Ts.npy')
Ts_flow_5L = np.load('Pilot_ox_base Flow 0.5\\Ts.npy')
Ts_flow_75 = np.load('Pilot_ox_base Flow 0.75\\Ts.npy')
Ts_flow_5H = np.load('Pilot_ox_base Flow 1.5\\Ts.npy')
Ts_flow_25 = np.load('Pilot_ox_base Flow 1.25\\Ts.npy')


# baseline plot
fig, ax = plt.subplots(ncols=2)
ax[0].plot(ts,Ts_base[-1,:],'tab:orange')
ax[0].set_xlabel('Time (s)')
ax[0].set_ylabel('Temperature (K)')
ax[0].set_title('Outlet Temperature During Oxidation Stage')
ax[1].plot(zs,Ts_base[3:-1,-1],'tab:orange')
ax[1].set_xlabel('Reactor Length (m)')
ax[1].set_title('Bed Temperature at End of Oxidation Stage')
#ax[1].yaxis.set_tick_params(labelbottom=True)
ax[0].text(-0.2,-0.1,'a)',transform=ax[0].transAxes,fontsize=14)
ax[1].text(-0.1,-0.1,'b)',transform=ax[1].transAxes,fontsize=14)


fig1, ax1 = plt.subplots(ncols=2)
ax1[0].plot(ts,Ts_flow_5L[-1,:],'tab:green',linestyle='dotted',label='0.113')#' kg/m$^2$/s')
ax1[0].plot(ts,Ts_flow_75[-1,:],'tab:orange',linestyle='dashed',label='0.170')# kg/m$^2$/s')
ax1[0].plot(ts,Ts_base[-1,:],'tab:blue',linestyle='solid',label='0.227')# kg/m$^2$/s')
ax1[0].plot(ts,Ts_flow_25[-1,:],'tab:red',linestyle='dashdot',label='0.284')# kg/m$^2$/s')
ax1[0].plot(ts,Ts_flow_5H[-1,:],'tab:purple',linestyle=(0, (3, 5, 1, 5, 1, 5)),label='0.341')# kg/m$^2$/s')
ax1[0].set_xlabel('Time (s)')
ax1[0].set_ylabel('Outlet Temperature (K)')
ax1[0].legend(title='Inlet Flux (kg/m$^2$/s)',title_fontsize=11,handlelength=2.5)

ax1[1].plot(zs,Ts_flow_5L[3:-1,-1],'tab:green',linestyle='dotted',label='0.113')# kg/m$^2$/s')
ax1[1].plot(zs,Ts_flow_75[3:-1,-1],'tab:orange',linestyle='dashed',label='0.170')#kg/m$^2$/s')
ax1[1].plot(zs,Ts_base[3:-1,-1],'tab:blue',linestyle='solid',label='0.227')# kg/m$^2$/s')
ax1[1].plot(zs,Ts_flow_25[3:-1,-1],'tab:red',linestyle='dashdot',label='0.284')# kg/m$^2$/s')
ax1[1].plot(zs,Ts_flow_5H[3:-1,-1],'tab:purple',linestyle=(0, (3, 5, 1, 5, 1, 5)),label='0.341')# kg/m$^2$/s')
ax1[1].set_xlabel('Reactor Length (m)')
#ax1[1].yaxis.set_tick_params(labelbottom=True)
ax1[1].set_ylabel('Final Temperature (K)')
#ax1[1].legend(title='Inlet mass flux',title_fontsize=11)
ax1[0].text(-0.2,-0.1,'a)',transform=ax1[0].transAxes,fontsize=14)
ax1[1].text(-0.1,-0.1,'b)',transform=ax1[1].transAxes,fontsize=14)


Ts_OC_3 = np.load('Pilot_ox_base OC 0.3\\Ts.npy')
Ts_OC_4 = np.load('Pilot_ox_base OC 0.4\\Ts.npy')
Ts_OC_6 = np.load('Pilot_ox_base OC 0.6\\Ts.npy')
Ts_OC_7 = np.load('Pilot_ox_base OC 0.7\\Ts.npy')

fig2, ax2 = plt.subplots(ncols=2)
ax2[0].plot(ts,Ts_OC_3[-1,:],'tab:green',linestyle='dotted',label='0.3')# kg Fe$_2$O$_3$/kg Red Mud')
ax2[0].plot(ts,Ts_OC_4[-1,:],'tab:orange',linestyle='dashed',label='0.4')# kg Fe$_2$O$_3$/kg Red Mud')
ax2[0].plot(ts,Ts_base[-1,:],'tab:blue',linestyle='solid',label='0.5')# kg Fe$_2$O$_3$/kg Red Mud')
ax2[0].plot(ts,Ts_OC_6[-1,:],'tab:red',linestyle='dashdot',label='0.6')# kg Fe$_2$O$_3$/kg Red Mud')
ax2[0].plot(ts,Ts_OC_7[-1,:],'tab:purple',linestyle=(0, (3, 5, 1, 5, 1, 5)),label='0.7')# kg Fe$_2$O$_3$/kg Red Mud')
ax2[0].set_xlabel('Time (s)')
ax2[0].set_ylabel('Outlet Temperature (K)')
ax2[0].legend(title='  OC Concentration\n(kg Fe$_2$O$_3$/kg Red Mud)',title_fontsize=11,handlelength=2.5)
#ax2[1].yaxis.set_tick_params(labelbottom=True)
ax2[0].text(-0.2,-0.1,'a)',transform=ax2[0].transAxes,fontsize=14)
ax2[1].text(-0.1,-0.1,'b)',transform=ax2[1].transAxes,fontsize=14)

ax2[1].plot(zs,Ts_OC_3[3:-1,-1],'tab:green',linestyle='dotted',label='0.3')# kg Fe$_2$O$_3$/kg Red Mud')
ax2[1].plot(zs,Ts_OC_4[3:-1,-1],'tab:orange',linestyle='dashed',label='0.4')# kg Fe$_2$O$_3$/kg Red Mud')
ax2[1].plot(zs,Ts_base[3:-1,-1],'tab:blue',linestyle='solid',label='0.5')# kg Fe$_2$O$_3$/kg Red Mud')
ax2[1].plot(zs,Ts_OC_6[3:-1,-1],'tab:red',linestyle='dashdot',label='0.6')# kg Fe$_2$O$_3$/kg Red Mud')
ax2[1].plot(zs,Ts_OC_7[3:-1,-1],'tab:purple',linestyle=(0, (3, 5, 1, 5, 1, 5)),label='0.7')# kg Fe$_2$O$_3$/kg Red Mud')
ax2[1].set_xlabel('Reactor Length (m)')
ax2[1].set_ylabel('Final Temperature (K)')
#plt.legend(title='OC Concentration',title_fontsize=11)


Ts_O2_15 =  np.load('Pilot_ox_base O2 0.15\\Ts.npy')
Ts_O2_165 = np.load('Pilot_ox_base O2 0.165\\Ts.npy')
Ts_O2_195 = np.load('Pilot_ox_base O2 0.195\\Ts.npy')
Ts_O2_21 = np.load('Pilot_ox_base O2 0.21\\Ts.npy')

fig3, ax3 = plt.subplots(ncols=2)
ax3[0].plot(ts,Ts_O2_15[-1,:],'tab:green',linestyle='dotted',label='0.15')
ax3[0].plot(ts,Ts_O2_165[-1,:],'tab:orange',linestyle='dashed',label='0.165')
ax3[0].plot(ts,Ts_base[-1,:],'tab:blue',linestyle='solid',label='0.18')
ax3[0].plot(ts,Ts_O2_195[-1,:],'tab:red',linestyle='dashdot',label='0.195')
ax3[0].plot(ts,Ts_O2_21[-1,:],'tab:purple',linestyle=(0, (3, 5, 1, 5, 1, 5)),label='0.21')
ax3[0].set_xlabel('Time (s)')
ax3[0].set_ylabel('Outlet Temperature (K)')
ax3[0].legend(title='O$_2$ Concentration',title_fontsize=11,handlelength=2.5)


ax3[1].plot(zs,Ts_O2_15[3:-1,-1],'tab:green',linestyle='dotted',label='0.15')
ax3[1].plot(zs,Ts_O2_165[3:-1,-1],'tab:orange',linestyle='dashed',label='0.165')
ax3[1].plot(zs,Ts_base[3:-1,-1],'tab:blue',linestyle='solid',label='0.18')
ax3[1].plot(zs,Ts_O2_195[3:-1,-1],'tab:red',linestyle='dashdot',label='0.195')
ax3[1].plot(zs,Ts_O2_21[3:-1,-1],'tab:purple',linestyle=(0, (3, 5, 1, 5, 1, 5)),label='0.21')
ax3[1].set_xlabel('Reactor Length (m)')
ax3[1].set_ylabel('Final Temperature (K)')
#ax3[1].yaxis.set_tick_params(labelbottom=True)
#ax3[1].legend(title='O$_2$ Concentration',title_fontsize=11)
ax3[0].text(-0.2,-0.1,'a)',transform=ax3[0].transAxes,fontsize=14)
ax3[1].text(-0.1,-0.1,'b)',transform=ax3[1].transAxes,fontsize=14)

# creating the tornado plot
from matplotlib.lines import Line2D
import pandas as pd
from highlight_text import ax_text

labels = np.char.array([
    'Inlet Mass Flux',#'\n -50% - +50%',
    'OC Concentration',#'\n 0.3 - 0.7 kg Fe$_2$O$_3$/kg ox Mud',
    'O$_2$ Inlet Fraction'])#'\n 0.15-0.21'])
# midpoint = 792.118 # MJ

# low_values = np.array([719.134,486.828,773.536]) # MJ

# high_values = np.array([796.386,1097.068,799.258]) # MJ

midpoint = 792.118/3.6 # kWh

low_values = np.array([719.134/3.6,486.828/3.6,773.536/3.6]) # kWh

high_values = np.array([796.386/3.6,1097.068/3.6,799.258/3.6]) # kWh

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
    plt.tight_layout()
    
    for y, low_value, high_value in zip(ys, low_values, high_values):
        low_width = midpoint - low_value
        high_width = high_value - midpoint
        
        plt.broken_barh([(low_value,low_width),(midpoint,high_width)], (y-0.35,0.7), facecolors=[color_low,color_high],
                        edgecolors=['black','black'], linewidth=0.5)
        offset = 25
        
        if high_value > low_value:
            x_high = midpoint + high_width + offset
            x_low = midpoint - low_width - offset
        else:
            x_high = midpoint + high_width - offset
            x_low = midpoint - low_width + offset
        #plt.text(x_high, y, str(round(high_value)), va='center', ha='center')
      #  plt.text(x_low, y, str(round(low_value)), va='center', ha = 'center')
        
    # for y, sub in zip(ys,['0.15-0.21','0.113-0.341 kg/m$^2$/s','0.3-0.7 kg Fe$_2$O$_3$/kg Red Mud']):
    #     plt.text( min(low_values)-100,y-0.22,sub,fontsize=10,ha='center',va='center',color='gray')
    
    plt.text(300/3.6,0,'O$_2$ Inlet Fraction', ha='center',va='center',color='black',fontsize=12)
    plt.text(300/3.6,2,'OC Concentration', ha='center', va='center',color='black',fontsize=12)
    plt.text(300/3.6,1,'Inlet Mass Flux',ha='center',va='center',color='black',fontsize=12)
    plt.axvline(midpoint,color='black',linewidth=1)
    ax_text(300/3.6,1-0.2, '<0.113> - <0.341> kg/m²/s',ha='center',va='center',color='black',fontsize=10,
            highlight_textprops=[{'color':'tab:blue','alpha':0.75},{'color':'tab:red','alpha':0.75}])
    ax_text(300/3.6,-0.2, '<0.15> - <0.21>',ha='center',va='center',color='black',fontsize=10,
            highlight_textprops=[{'color':'tab:blue','alpha':0.75},{'color':'tab:red','alpha':0.75}])
    
    ax_text(300/3.6,2-0.2, '<0.3> - <0.7> kg Fe₂O₃/kg Red Mud',ha='center',va='center',color='black',fontsize=10,
            highlight_textprops=[{'color':'tab:blue','alpha':0.75},{'color':'tab:red','alpha':0.75}])
    #plt.tight_layout()
    ax = plt.gca()
    ax.spines[['right','left','top']].set_visible(False)
    ax.set_yticks([])
    
    ax_text(x = midpoint, y = len(labels), s = title, color='black',fontsize=15,va='center',ha='center',
            highlight_textprops=[{'color':color_low,'fontweight':'bold'},{'color':color_high,'fontweight':'bold'}],ax=ax)
    ax.spines['left'].set_position(('data', min(low_values) - 25))
    plt.xlabel('Energy (kWh)')
   # plt.yticks(ys,labels,fontsize=14,ha='right')
    plt.xlim(450/3.6,1150/3.6)
    plt.ylim(-0.5,len(labels)-0.5)
    plt.tick_params(left=False)
    plt.text(792.118/3.6,2.65,'Baseline = 220',fontsize=10,ha='center',va='center')
    
    
    plt.show()
    return

tornado_chart(data['Labels'],midpoint,data['Low values'],data['High values'])


