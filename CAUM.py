
"""
Created on Tue Aug 21 22:24:29 2024
CAUM is an open-source software with a Python-based program that is developed to calculate
Chemical Ages of Uranium Minerals using electron probe microanalysis (EPMA) data.
The program consists of three main modules: "Single population", "Age extrapolation" and "Multi-population",
which are implemented by three sub-interfaces that can be invoked from three entries in the main user interface.
version 1.0
@author: Hao Song
@email: songhao@163.com
"""

import os
import math
import ctypes
import tkinter
import itertools
import numpy as np
import pandas as pd
import easygui as g
from tkinter import *
import seaborn as sns
from tkinter.ttk import *
import tkinter.font as tf
import ttkbootstrap as ttk
from ttkbootstrap import Style
import matplotlib.pyplot as plt
import tkinter.filedialog as fd
from scipy.stats import pearsonr
from tkinter.messagebox import *
from sklearn import linear_model
from scipy import stats, optimize
from Combopicker import Combopicker
from matplotlib.figure import Figure
from ttkbootstrap.constants import *
from scipy.spatial.distance import cdist
from scipy.stats.distributions import chi2
from tkinter.filedialog import asksaveasfilename
from Combopicker1 import Combopicker as Combopicker1
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
def Open():
    """ This module is used to obtain the data of the external excel
    table and display the obtained EPMA data path to the GUI interface."""
    global file
    file = g.fileopenbox(default="*.xlsx", filetypes=["*.xlsx"])
    path.insert(INSERT, file)
def Run():
    """ The module calculates the apparent age and error of each electron probe
     analysis point by using the content of UO2, ThO2 and PbO in the uranium
     mineral electron probe chemical composition analysis table input by the user
     and the corresponding error. At the same time, the module calculates the P-value
     and R-value between each oxide content. In addition, the module uses the calculated
     apparent age to draw the probability density curve of the apparent age."""
    def isFloat(x):
        """ The function is used to determine whether the error values of UO2, ThO2 and
        PbO contents input by the user are correct."""
        try:
            float(x)
            return True
        except:
            return False
    for i in tree1.get_children():# Clear list data in the ' The chemical Age ' module
        tree1.delete(i)
    for i in tree2.get_children():# Clear list data in the ' Coefficient of correlation ' module
        tree2.delete(i)
    for widget in frame3.winfo_children():# Clear the age kernel density curve in the ' The probability density plot ' module
        widget.destroy()
    if error_UO2.get().replace('.', '', 1).isdigit() and error_ThO2.get().replace('.', '', 1).isdigit() and error_PbO.get().replace('.', '', 1).isdigit():# To determine whether the error format of UO2, ThO2 and PbO entered by the user is correct.
        global df
        df = pd.read_excel(file)
        if 'Spot number' in df.columns:
            df = df.drop(columns=['Spot number'])
        else:
            pass
        # 检查每一列是否全部为数值
        non_numeric_columns = []
        for column in df.columns:
            if not pd.to_numeric(df[column], errors='coerce').notna().all():
                non_numeric_columns.append(column)
        if non_numeric_columns:#If there are non-numeric columns, use the pop-up window to remind the user
            mess = "Please check and modify the following non-numeric data：\n" + "\n".join(non_numeric_columns)
            tkinter.messagebox.showwarning("Non-numerical data warning", mess)
        else:
            if 'PbO' in df.columns and ('UO2' in df.columns or 'ThO2' in df.columns):
                if 0 in df['PbO'].values.tolist():
                    """ The content of PbO in the electron probe composition analysis table 
                    imported by the user cannot be zero. If the content of PbO is 0 or below 
                    the detection limit, the chemical age of the point cannot be calculated. 
                    The user should manually delete the analysis point data."""
                    tkinter.messagebox.showinfo('Message',
                                                'The percentage of PbO content at each analysis point can\'t be zero',
                                                parent=win)
                else:
                    if 'ThO2' in df:
                        pass
                    else:
                        example_list = [0, ] * len(df)
                        df['ThO2'] = np.array(example_list)  # ThO2 content inserted in dataframe

                    """ According to the content of UO2, ThO2 and PbO, the apparent age of each electron 
                    probe analysis point is calculated by iterative method."""
                    W_Th = 264  # The molecular weight of ThO2
                    W_U = 270  # The molecular weight of UO2
                    W_Pb = 222  # The molecular weight of PbO
                    a1 = 0.00000000098485  # Decay constant of U235
                    b1 = 0.000000000155125  # Decay constant of U238
                    r232 = 0.000000000049475  # Decay constant of Th232
                    ThO2 = df['ThO2']
                    UO2 = df['UO2']
                    PbO = df['PbO']
                    list_ThO2 = ThO2.values.tolist()  # series to list
                    list_UO2 = UO2.values.tolist()  # series to list
                    list_PbO = PbO.values.tolist()  # series to list
                    list_vaild = list(
                        zip(list_ThO2, list_UO2, list_PbO))  # The corresponding oxide package forms a list
                    list_age = []  # Create an empty list to store the age after iteration
                    for i in list_vaild:
                        i = list(i)
                        t = i[2] * 100 * 1e6
                        while ((i[0] / W_Th) * (math.exp(r232 * t) - 1) + (i[1] / W_U) * (
                                ((math.exp((a1) * t) + 137.88 * math.exp((b1) * t)) / 138.88) - 1)) > i[2] / W_Pb:
                            t = t - 1e4  # Iteration is repeated with 0.0001 Ma as a unit.
                        t = t * 1e-6  # The single-point chemical age results were converted to Ma as a unit
                        list_age.append(t)  # Put each single point chemical age into an empty list.
                    df['Age(Ma)'] = np.array(list_age)  # Inserting apparent age in dataframe
                    df['Age(Ma)'] = df['Age(Ma)'].apply(
                        lambda x: round(x, 2))  # The apparent age retains two decimal places

                    """ Based on the error transfer formula, the error of each apparent age is 
                    calculated according to the error of A, B and C input by the user."""
                    list_error_t = []  # 创建空列表存入迭代后年龄误差
                    if isFloat(error_UO2.get()) == True and isFloat(error_ThO2.get()) == True and isFloat(
                            error_PbO.get()) == True:
                        for i in list(zip(list_ThO2, list_UO2, list_PbO, list_age)):
                            # Calculation of UO2 transfer coefficient
                            K_UO2 = ((138.88 - math.exp(a1 * i[3]) - 137.88 * math.exp(b1 * i[3])) / (270 * 138.88)) / (
                                    (r232 * i[0] * math.exp(r232 * i[3]) / 264) + (
                                    a1 * i[1] * math.exp(a1 * i[3]) + 137.88 * b1 * i[1] * math.exp(b1 * i[3])) / (
                                            270 * 138.88))

                            # Calculation of ThO2 transfer coefficient
                            K_ThO2 = ((1 - math.exp(r232 * i[3])) / 264) / (
                                    (r232 * i[0] * math.exp(r232 * i[3]) / 264) + (
                                    a1 * i[1] * math.exp(a1 * i[3]) + 137.88 * b1 * i[1] * math.exp(b1 * i[3])) / (
                                            270 * 138.88))

                            # Calculation of PbO transfer coefficient
                            K_PbO = (1 / 222) / ((r232 * i[0] * math.exp(r232 * i[3]) / 264) + (
                                    a1 * i[1] * math.exp(a1 * i[3]) + 137.88 * b1 * i[1] * math.exp(b1 * i[3])) / (
                                                         270 * 138.88))

                            dev_UO2 = float(error_UO2.get()) / 100  # Get the value of dev _ UO2
                            dev_ThO2 = float(error_ThO2.get()) / 100  # Get the value of dev _ ThO2
                            dev_PbO = float(error_PbO.get()) / 100  # Get the value of dev _ PbO

                            # The absolute error of age Δt is calculated.
                            abs_t = K_UO2 * abs(i[1] * dev_UO2) + K_ThO2 * abs(i[0] * dev_ThO2) + K_PbO * abs(
                                i[2] * dev_PbO)
                            abs_t = abs_t * 1e-6  # Converted to Ma as a unit
                            list_error_t.append(abs_t)

                    df['Error(Ma)'] = np.array(list_error_t)
                    df['Error(Ma)'] = df['Error(Ma)'].apply(
                        lambda x: round(x, 2))  # The apparent age error retains two decimal places
                    df_datashow = df
                    df_datashow['ID'] = np.array(list(range(1, len(df_datashow) + 1)))
                    df = df.drop(['ID'], axis=1)
                    df_datashow['ID'] = df_datashow['ID'].astype(str)
                    data_show = df_datashow[['ID', 'UO2', 'ThO2', 'PbO', 'Age(Ma)', 'Error(Ma)']]
                    for i in np.array(data_show).tolist():
                        tree1.insert('', index=END, values=(i[0], i[1], i[2], i[3], i[4], i[5]))

                    """ Calculate the P-value and R-value between valid columns"""
                    df_datashow = df_datashow.drop(['ID'], axis=1)
                    if 'Total' in df_datashow:
                        df_datashow = df_datashow.drop(['Total'], axis=1)
                    if 'Age(Ma)' in df_datashow:
                        df_datashow = df_datashow.drop(['Age(Ma)'], axis=1)
                    if 'Error(Ma)' in df_datashow:
                        df_datashow = df_datashow.drop(['Error(Ma)'], axis=1)
                    correlation = df_datashow.corr()
                    result_P_R = []
                    for col1 in correlation.columns:
                        for col2 in correlation.columns:
                            if col1 != col2:
                                if df_datashow[col1].nunique() > 1 and df_datashow[col2].nunique() > 1:
                                    r, p = pearsonr(df_datashow[col1], df_datashow[col2])
                                    r = format(r, '.9f')  # The r-value retains 9 valid decimals.
                                    p = format(p, '.9f')  # The p-value retains 9 valid decimals.
                                    result_P_R.append((col1, col2, r, p))
                    sorted_data = sorted(result_P_R, key=lambda x: x[0] == 'PbO', reverse=True)
                    tag_colors = {}
                    for i in sorted_data:
                        element = i[0]
                        if element not in tag_colors:
                            if len(tag_colors) % 2 == 0:
                                tag_colors[element] = 'black'
                            else:
                                tag_colors[element] = 'black'
                        tree2.tag_configure(element, foreground=tag_colors[element])
                        tree2.insert('', index=END, text=i[0], values=(i[0], i[1], i[3], i[2]), tags=element)

                    """ Draw the distribution histogram of apparent age"""
                    # Create Matplotlib graphics
                    fig1 = Figure(figsize=(6.15, 5.05))
                    ax1 = fig1.add_subplot(111, facecolor='#DCDCDC')
                    ax1.grid(axis='x', linestyle='-', color='gray', linewidth=0.5)
                    # Draw a histogram
                    ax1.hist(list_age, density=False, color='cyan', edgecolor='gray')
                    # Add titles and labels
                    ax1.set_xlabel('Age(Ma)', fontdict=dict(weight='medium', family='Arial'))
                    ax1.set_ylabel('Number', fontdict=dict(weight='medium', family='Arial'))
                    ax2 = ax1.twinx()
                    # Add probability density curve
                    sns.kdeplot(list_age, color='red', ax=ax2, linewidth=3)
                    ax2.set_ylabel('')
                    ax2.set_yticks([])
                    ax2.spines['right'].set_visible(False)
                    for label in ax1.get_xticklabels():
                        label.set_fontname('Arial')
                    for label in ax1.get_yticklabels():
                        label.set_fontname('Arial')
                    # Embedding Matplotlib graphics into Tkinter
                    canvas1 = FigureCanvasTkAgg(fig1, master=frame3)
                    canvas1.draw()
                    canvas1.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=tkinter.YES)
                    toolbar1 = NavigationToolbar2Tk(canvas1, frame3)
                    toolbar1.update()
                    canvas1._tkcanvas.pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=tkinter.YES)
            else:
                tkinter.messagebox.showinfo('Message',
                                            '''The original file must contain columns named 'UO2(ThO2)' and 'PbO' ''',
                                            parent=win)
    else:
        tkinter.messagebox.showinfo("Message", "Please input the reasonable error value(UO2、ThO2 and PbO)")
def Reset():
    """ This module is used to clear all the content of the main
    interface to facilitate the next data processing"""
    path.delete('1.0', 'end')
    error_UO2.delete(0, END)
    error_ThO2.delete(0, END)
    error_PbO.delete(0, END)
    for i in tree1.get_children():
        tree1.delete(i)
    for i in tree2.get_children():
        tree2.delete(i)
    for widget in frame3.winfo_children():
        widget.destroy()
def Save_excel():
    """ This module is used to save the calculation results of apparent age and error."""
    try:
        df
        df.insert(0, 'Spot number', range(1, len(df) + 1))
        savefile = asksaveasfilename(filetypes=(("Excel files", "*.xlsx"), ("All files", "*.*")))
        df.to_excel(savefile + ".xlsx", index=False)
    except NameError:
        tkinter.messagebox.showinfo('Message', 'Please input the initial data', parent=win)

def weighted_age():
    """ This module implements the module that is designed to
    calculate the weighted mean chemical age of a single geological
    event."""
    def secection():
        """ This custom function is used to obtain external raw data(.xlsx)."""
        global file1
        res = g.fileopenbox(default="*.xlsx", filetypes=["*.xlsx"])
        file1 = res
        show_file.insert(INSERT, file1)
    def reset_input():
        """ This custom function is used to reset this sub-interface."""
        Result_age.delete(0, END)
        Result_error.delete(0, END)
        Result_MSWD.delete(0, END)
        Result_P.delete(0, END)
        show_file.delete(0, END)
        for widget in Weight_average.winfo_children():
            widget.destroy()
    def average_age():
        """ This custom function is used to calculate the the weighted average
        age, the error(σ), the MSWD (mean square of weighted deviation), and
        the P (probability)."""
        global df_1
        df_1 = pd.read_excel(file1)
        if 'Age(Ma)' in df_1.columns and 'Error(Ma)' in df_1.columns:
            if valii.get()==2:  # When the error type is perecent(%)
                df_1['Error(Ma)']=df_1.apply(lambda row: row['Age(Ma)'] * row['Error(Ma)'], axis=1)
            else:               # When the error type is abs.(ppm)
                pass
            list_age=df_1['Age(Ma)'].tolist()
            list_error_t = df_1['Error(Ma)'].tolist()
            age_up = []
            age_down = []
            total_error = []
            mswd = []
            for i in list(zip(list_age, list_error_t)):
                if vali.get() == 1: # When the error type is 1 sigma
                    age_up.append(i[0] / (i[1] ** 2))
                    age_down.append(1 / (i[1] ** 2))
                    total_error.append(1 / (i[1] ** 2))
                    mswd.append(
                        ((i[0] - sum(list_age) / len(list_age)) ** 2 / (i[1] ** 2)) / (len(list_age) - 1))
                else: # When the error type is 2 sigma
                    age_up.append(i[0] / ((i[1] / 2) ** 2))
                    age_down.append(1 / ((i[1] / 2) ** 2))
                    total_error.append(1 / ((i[1] / 2) ** 2))
                    mswd.append(
                        ((i[0] - sum(list_age) / len(list_age)) ** 2 / ((i[1] / 2) ** 2)) / (len(list_age) - 1))
            average_t = sum(age_up) / sum(age_down)
            total_error = 2 / (sum(total_error) ** 0.5)
            MSWD = sum(mswd)
            average_t = '{:.2f}'.format(average_t)      # Calcute the weighted average age (T)
            total_error = '{:.2f}'.format(total_error)  # Calcute the standard deviation (σ)
            MSWD = '{:.2f}'.format(MSWD)                # Calcute the  mean square of weighted deviation (MSWD)
            polit = chi2.sf(float(MSWD) * (len(list_age) - 1), len(list_age) - 1) # Calcute the P (probability)
            Result_age.delete(0, END)
            Result_error.delete(0, END)
            Result_MSWD.delete(0, END)
            Result_P.delete(0, END)
            for widget in Weight_average.winfo_children():
                widget.destroy()
            Result_age.insert(INSERT, average_t)        # Show the weighted average age (T)
            Result_error.insert(INSERT, total_error)    # Show the standard deviation (σ)
            Result_MSWD.insert(INSERT, MSWD)            # Show the  mean square of weighted deviation (MSWD)
            Result_P.insert(INSERT, polit)              # Show the P (probability)

            """ Using the calculated apparent age and error, the weighted average 
            age map is drawn and displayed in the “Weighted Average” module."""
            f = Figure(figsize=(6, 5), dpi=110)
            a = f.add_subplot(111)
            y = df_1['Age(Ma)'].values.tolist()
            x = list(range(1, len(y) + 1))
            yerr = df_1['Error(Ma)'].values.tolist()
            a.errorbar(x, y, yerr, fmt='co', mec='red', elinewidth=8, ecolor='red', mfc='red', ms=1)
            a.axhline(y=float(average_t), c='green', lw=2, zorder=0)  # raw the average age straight line
            a.tick_params(axis='x', bottom=False, colors='white')
            a.tick_params(axis='y', tickdir='in')  # Y轴刻度指向图内
            a.set_ylabel('Age(Ma)', fontdict=dict(weight='medium', family='Arial'))  # Set the Y-axis label
            a.set_xlim(0, len(y) + 2)  # 设置X轴范围
            a.set_ylim(min(y) - 2 * yerr[y.index(min(y))], max(y) + 2 * yerr[y.index(max(y))])  # Set the Y-axis range
            a.set_facecolor('#DCDCDC')  # Set the background color of the weighted graph
            a.grid(axis='y')

            # Set the upper right corner annotation ( weighted graph drawing type : 1σ or 2σ )
            if vali.get() == 1:
                a.set_title('box heights are ' + str(1) + 'σ', loc='right',
                            fontdict=dict(weight='medium', family='Arial', size=10), pad=5)
            else:
                a.set_title('box heights are ' + str(2) + 'σ', loc='right',
                            fontdict=dict(weight='medium', family='Arial', size=10), pad=5)
            # Insert age calculation results in the weighted average age map
            a.text(0.83, 0.92, s='Mean=' + str(average_t) + '±' + str(total_error) + ' Ma' + '\nMSWD=' + str(
                MSWD), fontdict=dict(weight='medium', family='Arial'),
                   bbox={'facecolor': 'white', 'edgecolor': 'black', 'pad': 0.5, 'boxstyle': 'round',
                         'alpha': 0.8}, ha='center',
                   va='center', transform=a.transAxes)
            canvas = FigureCanvasTkAgg(f, master=Weight_average)
            canvas.draw()
            canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=tkinter.YES)
            toolbar = NavigationToolbar2Tk(canvas, Weight_average)
            toolbar.update()
            canvas._tkcanvas.pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=tkinter.YES)
            def on_key_event(event):
                key_press_handler(event, canvas, toolbar)
            canvas.mpl_connect('key_press_event', on_key_event)
        else:
            tkinter.messagebox.showinfo('Message', 'The initial table must contain columns named Age(Ma) and Error(Ma)', parent=t2)
    t2 = Toplevel()
    t2.title('Single population')
    t2.geometry(s_center)
    t2.geometry('988x845+350+220')
    t2.resizable(False, False)
    Operation = tkinter.LabelFrame(t2, text='Operation', width=303, height=840, font='Times 17 bold')
    Operation.grid(row=0, column=0, padx=5)
    Operation.grid_propagate(False)
    Weight_average = tkinter.LabelFrame(t2, text='Weighted Average', width=670, height=610, font='Times 17 bold')
    Weight_average.grid(row=0, column=1, sticky=N)
    Weight_average.grid_propagate(False)
    Data_input = tkinter.LabelFrame(Operation, text='Input', width=290, height=150, font='Times 12 bold')
    Data_input.grid(row=0, column=0, padx=0, pady=5, rowspan=6, columnspan=2, sticky=N)
    Data_input.grid_propagate(False)
    Result = tkinter.LabelFrame(Operation, text='Result', width=290, height=175, font='Times 12 bold')
    Result.grid(row=10, column=0, padx=5, pady=5, rowspan=6, columnspan=2, sticky=N)
    Result.grid_propagate(False)
    Type_error = tkinter.LabelFrame(Data_input, text='Type(error)', width=280, height=110, font='Times 12 bold')
    Type_error.grid(row=4, column=0, padx=5, pady=0, rowspan=6, columnspan=2, sticky=N)
    Type_error.grid_propagate(False)
    # Set the error type(1 sigma or 2 sigma)
    vali = IntVar()
    vali.set('2')  # The default error type is 2 sigma
    tkinter.Radiobutton(Type_error, variable=vali, value='1', text='1sigma', font='Times 15 bold').grid(row=0, column=0,padx=15, pady=5)
    tkinter.Radiobutton(Type_error, variable=vali, value='2', text='2sigma', font='Times 15 bold').grid(row=1, column=0,pady=5)
    # Set the error type(absolute or percent)
    valii = IntVar()
    valii.set('1')  # The abs is selected as the error type by default.
    tkinter.Radiobutton(Type_error, variable=valii, value='1', text='abs.(ppm)', font='Times 15 bold').grid(row=0,column=1,sticky=W)
    tkinter.Radiobutton(Type_error, variable=valii, value='2', text='percent(%)', font='Times 15 bold').grid(row=1,column=1,sticky=W)
    tkinter.Button(Operation, text='Open', width=12, font='Times 15 ', command=secection).grid(row=6, column=0,columnspan=2, padx=5,pady=5,sticky=E)
    show_file = tkinter.Entry(Operation, relief='groove', width=18)
    show_file.grid(row=6, column=0, pady=5, padx=5, columnspan=2, sticky=W + N + S)
    tkinter.Button(Operation, text='Run', width=26, font='Times 15 ', command=average_age).grid(row=7, column=0, pady=5,padx=5)
    tkinter.Button(Operation, text='Reset', width=26, font='Times 15 ', command=reset_input).grid(row=8, column=0,pady=5, padx=5)
    tkinter.Label(Result, text='Age(Ma):', font='Times 15 bold').grid(row=10, column=0, padx=5, pady=5, sticky=E)
    tkinter.Label(Result, text='Error(Ma):', font='Times 15 bold').grid(row=11, column=0, padx=5, pady=5, sticky=E)
    tkinter.Label(Result, text='MSWD:', font='Times 15 bold').grid(row=12, column=0, padx=5, pady=5, sticky=E)
    tkinter.Label(Result, text='P(Probability):', font='Times 15 bold').grid(row=13, column=0, padx=5, pady=5, sticky=W)
    Result_age = tkinter.Entry(Result, relief='groove', width=15)
    Result_age.grid(row=10, column=1, pady=5, padx=5, sticky=W)
    Result_error = tkinter.Entry(Result, relief='groove', width=15)
    Result_error.grid(row=11, column=1, pady=5, padx=5, sticky=W)
    Result_MSWD = tkinter.Entry(Result, relief='groove', width=15)
    Result_MSWD.grid(row=12, column=1, pady=5, padx=5, sticky=W)
    Result_P = tkinter.Entry(Result, relief='groove', width=15)
    Result_P.grid(row=13, column=1, pady=5, padx=5, sticky=W)
    tkinter.Button(Operation, text='Exit', width=26, font='Times 15', command=lambda: t2.destroy()).grid(row=16,column=0,pady=5, padx=5,sticky=W)
def estimate_age():
    """ This module is used to calculate the crystallization age of the uranium
    mineral through data extrapolation if there is a demonstrable linear relationship
    between PbO contents or chemical ages and impurities."""
    def secection():
        """ This module is used to import electron probe data and preprocess tables."""
        global file2
        file2 = g.fileopenbox(default="*.xlsx", filetypes=["*.xlsx"])
        show_file.insert(INSERT, file2)
        df_estimate = pd.read_excel(file2)
        if 'Spot number' in df_estimate.columns:
            df_estimate = df_estimate.drop(columns=['Spot number'])
        else:
            pass
        if 'Total' in df_estimate:
            df_estimate = df_estimate.drop(['Total'], axis=1)  # 去除名为‘Total’的列
        if 'Age(Ma)' in df_estimate:
            df_estimate = df_estimate.drop(['Age(Ma)'], axis=1)  # 去除名为‘Age(Ma)’的列
        if 'Error(Ma)' in df_estimate:
            df_estimate = df_estimate.drop(['Error(Ma)'], axis=1)  # 去除名为‘Error(Ma)’的列
        X_element.values = list(df_estimate.keys())
        Y_element.values = list(df_estimate.keys()) + ['Age(Ma)']
    def clear():
        """ This function is used to reset the "Age extrapolation" interface."""
        show_file.delete(0, END)  # Clears the display path text box
        plot_title.delete(0, END)  # Clear the graphic title
        X_plot_title.delete(0, END)  # Clear the X-axis label
        Y_plot_title.delete(0, END)  # Clear the Y-axis label
        X_element.delete(0, END)
        Y_element.delete(0, END)
        X_element.values = []
        Y_element.values = []
        set_bands.delete(0, END)
        result_slope.delete(0, END)
        result_equation.delete(0, END)
        result_intercept.delete(0, END)
        result_R2.delete(0, END)
        # Clear the data in the treeview
        for widget in labelframe_Age.winfo_children():
            widget.destroy()
    def plot():
        """ The module is used to calculate the apparent age of each electron probe
        analysis point. The correlation between oxides is used to calculate the
        crystallization age of uranium minerals and draw the correlation diagram."""
        for widget in labelframe_Age.winfo_children():
            widget.destroy()
        result_slope.delete(0, END)
        result_equation.delete(0, END)
        result_intercept.delete(0, END)
        result_R2.delete(0, END)
        if len(show_file.get()) == 0:
            tkinter.messagebox.showinfo('Message', 'Please selcet the file(.xlsx) !', parent=t3)
        else:
            df_plot = pd.read_excel(file2)
            if 'Spot number' in df_plot.columns:
                df_plot = df_plot.drop(columns=['Spot number'])
            else:
                pass
            example_list = [0, ] * len(df_plot)
            if 'ThO2' in df_plot:
                pass
            else:
                df_plot['ThO2'] = np.array(example_list)
            if 'Age(Ma)' in df_plot:
                pass
            else:
                # The apparent age of each point is calculated according to the iterative method.
                W_Th = 264
                W_U = 270
                W_Pb = 222
                ThO2 = df_plot['ThO2']
                UO2 = df_plot['UO2']
                PbO = df_plot['PbO']
                a1 = 0.00000000098485
                b1 = 0.000000000155125
                r232 = 0.000000000049475
                list_ThO2 = ThO2.values.tolist()
                list_UO2 = UO2.values.tolist()
                list_PbO = PbO.values.tolist()
                list_vaild = list(zip(list_ThO2, list_UO2, list_PbO))
                list_age = []
                for i in list_vaild:
                    i = list(i)
                    t = i[2] * 100 * 1e6
                    while ((i[0] / W_Th) * (math.exp(r232 * t) - 1) + (i[1] / W_U) * (
                            ((math.exp((a1) * t) + 137.88 * math.exp((b1) * t)) / 138.88) - 1)) > i[2] / W_Pb:
                        t = t - 1e4
                    t = t * 1e-6
                    list_age.append(t)
                df_plot['Age(Ma)'] = np.array(list_age)
            if len(X_element.get()) == 0 or len(Y_element.get()) == 0:
                tkinter.messagebox.showinfo('Message', 'x(y) parameter cannot be null', parent=t3)
            else:
                X_list = X_element.get().split('+')  # Get the list of oxides that need to be drawn for the X axis
                Y_list = Y_element.get().split('+')  # Get the list of oxides that need to be drawn for the Y axis
                x = sum([df_plot[i] for i in X_list]).values
                y = sum([df_plot[i] for i in Y_list]).values
                # Linefit
                b, a = np.polyfit(x, y, 1)
                x_pred = np.linspace(min(x), max(x), 100)
                y_pred = b * x_pred + a
                # Draw scatter plots and fit straight lines
                fig, ax = plt.subplots()
                ax.plot(x, y, marker=str(shape_spot.get()), color=str(color_spot.get()), linestyle='',
                        label='The measured value')
                ax.plot(x_pred, y_pred, 'k', label='Linear fitting')
                # 计算P值与R2
                P_value = (stats.pearsonr(x, y))[1]
                r2_value = ((stats.pearsonr(x, y))[0]) ** 2
                ax.set_xlabel(X_plot_title.get(), fontdict=dict(weight='medium', family='Arial'))
                ax.set_ylabel(Y_plot_title.get(), fontdict=dict(weight='medium', family='Arial'))
                ax.tick_params(axis='x', tickdir='in')  # X-axis scale pointer inward
                ax.tick_params(axis='y', tickdir='in')  # Y-axis scale pointer inward
                ax.set_title(plot_title.get(), fontdict=dict(weight='medium', family='Arial'),
                             size=10, loc='center')
                ax.set_title('Slope: {:.7f}'.format(b) + '\nIntercept: {:.7f}'.format(a),
                             fontdict=dict(weight='medium', family='Arial', size=10), loc='left')
                ax.set_title('R2: {:.7f}'.format(r2_value) + '\nP: {:.7f}'.format(P_value),
                             fontdict=dict(weight='medium', family='Arial', size=10), loc='right')
                result_slope.insert(INSERT, b)
                result_intercept.insert(INSERT, a)
                result_R2.insert(INSERT, r2_value)
                if a > 0:
                    result_equation.insert(INSERT, 'y=' + str(round(b, 3)) + 'x+' + str(round(a, 3)))
                else:
                    result_equation.insert(INSERT, 'y=' + str(round(b, 3)) + 'x' + str(round(a, 3)))
                # Modify the font of the X-axis scale to Times New Roman
                for label in ax.get_xticklabels():
                    label.set_fontname('Arial')
                # Modify the Y-axis scale font to Arial
                for label in ax.get_yticklabels():
                    label.set_fontname('Arial')
                if Confidence_Bands.get():
                    bands_params = set_bands.get()
                    try:
                        if float(bands_params) >= 0 and float(bands_params) <= 1:
                            alpha = 1 - float(bands_params)
                            limit = (1 - alpha) * 100
                            residuals = y - (b * x + a)
                            std_error = np.std(residuals)
                            critical_value = stats.t.isf(alpha / 2, len(x) - 2)

                            if Confidence_Bands.get():  # If the confidence interval is selected
                                # Calculate the boundary value of the confidence interval of the Y-axis intercept
                                a_std_error = std_error * np.sqrt(
                                    1 / len(x) + np.mean(x) ** 2 / np.sum((x - np.mean(x)) ** 2))
                                a_lower_bound = a - critical_value * a_std_error
                                a_upper_bound = a + critical_value * a_std_error
                                print("Lower bound of y-intercept:", a_lower_bound)
                                print("Upper bound of y-intercept:", a_upper_bound)
                                # Draw Confidence interval
                                se_fit = std_error * np.sqrt(
                                    1 / len(x) + (x_pred - np.mean(x)) ** 2 / np.sum((x - np.mean(x)) ** 2))
                                ax.fill_between(x_pred, y_pred - critical_value * se_fit,
                                                y_pred + critical_value * se_fit,
                                                alpha=0.3, color='red',
                                                label='Confidence interval ({0:.1f}%)'.format(limit))
                                ax.set_title(
                                    'Slope: {:.7f}'.format(b) + '\nIntercept: {:.7f}'.format(a) + '±{:.2f}'.format(
                                        critical_value * a_std_error),
                                    fontdict=dict(weight='medium', family='Arial', size=10), loc='left')

                            ax.legend(loc='best',
                                      prop={'weight': 'medium', 'family': 'Arial', 'size': 'medium'})
                            # Create a Figure object and put it in the Tkinter window
                            fig_canvas = FigureCanvasTkAgg(fig, master=labelframe_Age)
                            fig_canvas.draw()
                            fig_canvas.get_tk_widget().pack()
                            # Create a toolbar object and put it in the Tkinter window
                            toolbar = NavigationToolbar2Tk(fig_canvas, labelframe_Age)
                            toolbar.update()
                            fig_canvas.get_tk_widget().pack()
                        else:
                            tkinter.messagebox.showinfo('Message',
                                                        'Confidence level must be a single numeric value greater than 0 and less than 1!',
                                                        parent=t3)
                    except ValueError:
                        tkinter.messagebox.showinfo('Message',
                                                    'Confidence level must be a single numeric value greater than 0 and less than 1!',
                                                    parent=t3)
                else:
                    ax.legend(loc='best', prop={'weight': 'medium', 'family': 'Arial', 'size': 'medium'})
                    # Create a Figure object and put it in the Tkinter window
                    fig_canvas = FigureCanvasTkAgg(fig, master=labelframe_Age)
                    fig_canvas.draw()
                    fig_canvas.get_tk_widget().pack()
                    # Create a toolbar object and put it in the Tkinter window
                    toolbar = NavigationToolbar2Tk(fig_canvas, labelframe_Age)
                    toolbar.update()
                    fig_canvas.get_tk_widget().pack()
    t3 = Toplevel()
    t3.title('Age extrapolation')
    t3.geometry(s_center)
    t3.geometry('910x720+370+240')
    t3.resizable(True, True)
    # Set the Operation tag box
    Operation = tkinter.LabelFrame(t3, text='Operation', width=250, height=700, font='Times 17 bold')
    Operation.grid(row=0, column=0, padx=5, sticky=N + S)
    Operation.grid_propagate(False)
    # Set the import file label box
    Data_input = tkinter.LabelFrame(Operation, text='Input', width=240, height=70, font='Times 15 bold')
    Data_input.grid(row=0, column=0, padx=5, pady=1, rowspan=2, columnspan=2, sticky=N + W)
    Data_input.grid_propagate(False)
    # Set the display file path text box
    show_file = tkinter.Entry(Data_input, relief='groove', width=15)
    show_file.grid(row=0, column=0, pady=5, padx=5, sticky=N + S)
    tkinter.Button(Data_input, text='Open', width=9, font='Times 15 ', command=secection).grid(row=0, column=1, padx=7,                                                                                  pady=5)
    # Set the sub-tag box of the drawing area
    Element = tkinter.LabelFrame(Operation, text='Selection', width=240, height=440, font='Times 15 bold')
    Element.grid(row=7, column=0, padx=5, pady=1, rowspan=6, columnspan=2, sticky=N)
    Element.grid_propagate(False)
    # Set the title of the fitting line
    tkinter.Label(Element, text='Title:', font='Times 15 bold').grid(row=0, column=0, pady=5, padx=5, sticky=N + S)
    plot_title = tkinter.Entry(Element, relief='groove', width=19)
    plot_title.grid(row=0, column=1, pady=5, padx=5, sticky=N + S + W)
    # Set the X-axis parameters
    tkinter.Label(Element, text='X:', font='Times 15 bold').grid(row=1, column=0, pady=5, padx=5, sticky=N + S)
    X_element = Combopicker(Element, entrywidth=18)
    X_element.grid(row=1, column=1, pady=5, padx=5, sticky=N + S + W)
    # Set the Y-axis parameters
    tkinter.Label(Element, text='Y:', font='Times 15 bold').grid(row=2, column=0, pady=5, padx=5, sticky=N + S)
    Y_element = Combopicker(Element, entrywidth=18)
    Y_element.grid(row=2, column=1, pady=5, padx=5, sticky=N + S + W)
    # Set the X-axis label
    tkinter.Label(Element, text='X Axes:', font='Times 15 bold').grid(row=3, column=0, pady=5, padx=5, sticky=N + S)
    X_plot_title = tkinter.Entry(Element, relief='groove', width=19)
    X_plot_title.grid(row=3, column=1, pady=5, padx=5, sticky=N + S + W)
    # Set the Y-axis label
    tkinter.Label(Element, text='Y Axes:', font='Times 15 bold').grid(row=4, column=0, pady=5, padx=5, sticky=N + S)
    Y_plot_title = tkinter.Entry(Element, relief='groove', width=19)
    Y_plot_title.grid(row=4, column=1, pady=5, padx=5, sticky=N + S + W)
    # Set the coordinate point color
    tkinter.Label(Element, text='Color:', font='Times 15 bold').grid(row=5, column=0, pady=10, padx=5, sticky=N + S)
    color_spot = Combobox(Element, state='readonly', width=16)
    color_spot['values'] = (
        'red', 'black', 'green', 'blue', 'yellow', 'purple', 'white', 'brown', 'orange', 'cyan', 'magenta', 'pink',
        'slategray', 'light', 'tan', 'grey')
    color_spot.current(0)
    color_spot.grid(row=5, column=1, pady=5, padx=5, sticky=W)
    # Set the coordinate point shape
    tkinter.Label(Element, text='Shape:', font='Times 15 bold').grid(row=6, column=0, pady=10, padx=5, sticky=N + S)
    shape_spot = Combobox(Element, state='readonly', width=16)
    shape_spot['values'] = ('*', 'o', '^', 'v', '<', '>', 's', 'P', 'p', 'H', 'h', 'D', 'd', 'X')
    shape_spot.current(0)
    shape_spot.grid(row=6, column=1, pady=5, padx=5, stick=W)
    # Set the confidence interval
    def change_bands_inputstate():
        if Confidence_Bands.get() or Prediction_Bands.get():
            set_bands.config(state=tkinter.NORMAL)
        else:
            set_bands.config(state=tkinter.DISABLED)
    Bands = tkinter.LabelFrame(Element, text='Evaluation of confidence', width=230, height=100, font='Times 15 bold')
    Bands.grid(row=7, column=0, padx=5, pady=1, rowspan=3, columnspan=2, sticky=N + W + E)
    Bands.grid_propagate(False)
    Confidence_Bands = tkinter.BooleanVar()
    Prediction_Bands = tkinter.BooleanVar()
    Confidence_Bands.set(False)
    Prediction_Bands.set(False)
    checkbox1_widget = tkinter.Checkbutton(Bands, text="     Confidence Bands", font='Times 15 bold', variable=Confidence_Bands, command=change_bands_inputstate)
    checkbox1_widget.grid(row=0, column=0, padx=10, sticky=N + S)
    tkinter.Label(Bands, text='Confidence Level:', font='Times 15 bold').grid(row=2, column=0, pady=5, padx=5,sticky=N + S + W)
    set_bands = tkinter.Entry(Bands, relief='groove', width=7, state=tkinter.DISABLED)
    set_bands.grid(row=2, column=0, pady=5, padx=5, sticky=N + S + E)
    # Set the image button
    tkinter.Button(Operation, text='Plot', width=9, font='Times 15 ', command=plot).grid(row=16, column=0, pady=5,padx=8, columnspan=2, sticky=W)
    # Settings Clear all displayed button
    tkinter.Button(Operation, text='Reset', width=9, font='Times 15 ', command=clear).grid(row=16, column=0, pady=5,padx=8, columnspan=2,sticky=E)
    # Set the exit button
    tkinter.Button(Operation, text='Exit', width=21, font='Times 15 ', command=lambda: t3.destroy()).grid(row=17,column=0,pady=5,padx=5,columnspan=2)
    # Set the age projection drawing area
    labelframe_Age = tkinter.LabelFrame(t3, text='Diagram of correlation', width=640, height=550, font='Times 17 bold')
    labelframe_Age.grid(row=0, column=1, padx=2, sticky=N + W)
    labelframe_Age.grid_propagate(False)
    # Set the result display area
    labelframe_eq = tkinter.LabelFrame(t3, text='Result', width=640, height=130, font='Times 17 bold')
    labelframe_eq.grid(row=0, column=1, padx=2, sticky=S + W)
    labelframe_eq.grid_propagate(False)
    tkinter.Label(labelframe_eq, text='Slope:', font='Times 17 bold').grid(row=0, column=0, pady=10, padx=5)
    result_slope = tkinter.Entry(labelframe_eq, relief='groove', width=26)
    result_slope.grid(row=0, column=1, pady=10, padx=5)
    tkinter.Label(labelframe_eq, text='Equation:', font='Times 17 bold').grid(row=1, column=0, pady=10, padx=5)
    result_equation = tkinter.Entry(labelframe_eq, relief='groove', width=26)
    result_equation.grid(row=1, column=1, pady=10, padx=5)
    tkinter.Label(labelframe_eq, text='Intercept:', font='Times 17 bold').grid(row=0, column=2, pady=10, padx=5)
    result_intercept = tkinter.Entry(labelframe_eq, relief='groove', width=26)
    result_intercept.grid(row=0, column=4, pady=10, padx=5)
    tkinter.Label(labelframe_eq, text='R2:', font='Times 17 bold').grid(row=1, column=2, pady=10, padx=5)
    result_R2 = tkinter.Entry(labelframe_eq, relief='groove', width=26)
    result_R2.grid(row=1, column=4, pady=10, padx=5)
def Mix_age():
    """ This module is used to calculate multiple age populations. """
    def count_sub_age_domain(data, n, sigma):
        """ Through the user input apparent age and error data set,
        age domain number and error type, the module can use the least
        squares extension model to calculate the best clustering model
        under different N values, and calculate the weighted average age
        of each sub-age domain(age,error,MSWD and P)."""
        data = sorted(data, key=lambda x: x[0])  # The data are arranged from small to large.
        # Generate index combination of partition scheme
        combinations = itertools.combinations(range(1, len(data)), n - 1)
        result = []#Used to store calculation results
        # Divided according to the index combination
        for comb in combinations:
            groups = []
            current_index = 0
            for i in comb:
                group = data[current_index:i]
                groups.append(group)
                current_index = i
            # The remaining part is taken as the last category.
            last_group = data[current_index:]
            groups.append(last_group)
            # Only the partition scheme containing non-empty list is added to the result.
            if all(groups):
                result.append(groups)
        def countssa(list1):  # Calculate each sub-age domain SSE
            age_up = []
            age_down = []
            for i in list1:
                age_up.append(i[0] / (i[1] ** 2))
                age_down.append(1 / (i[1] ** 2))
            av_t = sum(age_up) / sum(age_down)
            sub_s = []
            for i in list1:
                sub_s.append(((i[0] - av_t) / i[1]) ** 2)
            S = sum(sub_s)
            return S
        for i in result:
            for j in i:
                j.append(countssa(j))
        sub = []
        for i in result:
            sub_sse = 0
            for j in i:
                sub_sse += j[-1]  # Accumulate the last element
            sub.append(sub_sse)
        vaild_index = sub.index(min(sub))
        vaild_age_domain = result[vaild_index]  # Effectively distinguished sub-age results
        # Calculate the age of each sub-age domain
        total_result = []
        for i in vaild_age_domain:
            del i[-1]
            averageup = []
            averagedown = []
            suberror = []
            Mswd = []
            for j in i:
                if sigma == 1:
                    averageup.append(j[0] / (j[1] ** 2))
                    averagedown.append(1 / (j[1] ** 2))
                    suberror.append(1 / (j[1] ** 2))
                else:
                    averageup.append(j[0] / ((j[1] / 2) ** 2))
                    averagedown.append(1 / ((j[1] / 2) ** 2))
                    suberror.append(1 / ((j[1] / 2) ** 2))
            sub_average = sum(averageup) / sum(averagedown)
            sub_error = 2 / (sum(suberror) ** 0.5)
            if len(i) == 1:
                Mswd.append(0)
            else:
                for k in i:
                    if sigma == 1:
                        Mswd.append(((k[0] - sub_average) / k[1]) ** 2 / (len(i) - 1))
                    else:
                        Mswd.append(((k[0] - sub_average) / (k[1] / 2)) ** 2 / (len(i) - 1))
            MSWD = sum(Mswd)
            po =chi2.sf(float(MSWD) * (len(i) - 1), len(i) - 1)
            po = "{:.3e}".format(po)
            total_result.append((sub_average, sub_error,MSWD, len(i), po))
        totalS=[]
        for i in total_result:
            totalS.append(i[2]*i[3]-i[2])
        totalP = chi2.sf(sum(totalS), len(data)-n)
        total_result.append((sum(totalS), totalP))
        return total_result
    def secection():
        """ This module is used to import electron probe data."""
        global file3
        res = g.fileopenbox(default="*.xlsx", filetypes=["*.xlsx"])
        file3 = res
        show_file.insert(INSERT, file3)
    def clear_all():
        """ This module is used for reset."""
        show_file.delete(0, END)  # Clear the file path display box
        # Clear all Age results
        for widget in Result.winfo_children():
            widget.destroy()
        tkinter.Label(Result, text='------------------------Calculated------------------------', font='Times 25 bold').grid(row=0, column=0,padx=10,pady=10,columnspan=5,sticky=W + E)
        tkinter.Label(Result, text='Age(Ma)', font='Times 15 bold').grid(row=1, column=0, padx=10, pady=5)
        tkinter.Label(Result, text='Error(Ma)', font='Times 15 bold').grid(row=1, column=1, padx=10, pady=5)
        tkinter.Label(Result, text='MSWD', font='Times 15 bold').grid(row=1, column=2, padx=10, pady=5)
        tkinter.Label(Result, text='P', font='Times 15 bold').grid(row=1, column=3, padx=10, pady=5)
        tkinter.Label(Result, text='N', font='Times 15 bold').grid(row=1, column=4, padx=10, pady=5)
        for widget in SSE.winfo_children():
            widget.destroy()
        tkinter.Label(SSE, text='P_total:', font='Times 18 bold').grid(row=1, column=0, padx=10, rowspan=1,sticky=S + E)
    def count():
        df_mix = pd.read_excel(file3)  # The table must contain age and error.
        if 'Age(Ma)' and 'Error(Ma)' in df_mix:
            if type_content.get() == 2:
                df_mix['Error(Ma)'] = df_mix['Error(Ma)'] * df_mix['Age(Ma)']
            else:# The error type is absolute value.
                pass
            # Clear all age results
            for widget in Result.winfo_children():
                widget.destroy()
            # Add the title of the display box
            tkinter.Label(Result, text='------------------------Calculated------------------------',
                          font='Times 25 bold').grid(row=0, column=0,padx=10,pady=10,columnspan=5,sticky=W + E)
            tkinter.Label(Result, text='Age(Ma)', font='Times 15 bold').grid(row=1, column=0, padx=10, pady=5)
            tkinter.Label(Result, text='Error(Ma)', font='Times 15 bold').grid(row=1, column=1, padx=10, pady=5)
            tkinter.Label(Result, text='MSWD', font='Times 15 bold').grid(row=1, column=2, padx=10, pady=5)
            tkinter.Label(Result, text='P', font='Times 15 bold').grid(row=1, column=3, padx=10, pady=5)
            tkinter.Label(Result, text='N', font='Times 15 bold').grid(row=1, column=4, padx=10, pady=5)
            for widget in SSE.winfo_children():
                widget.destroy()
            tkinter.Label(SSE, text='P_total:', font='Times 18 bold').grid(row=1, column=0, padx=10, rowspan=1,sticky=S + E)
            input_data = df_mix[['Age(Ma)', 'Error(Ma)']].values.tolist()
            print_data = count_sub_age_domain(input_data, int(spin.get()), int(type_sigma.get()))
            tkinter.Label(SSE, text="{:.3e}".format(print_data[-1][1]), font='Times 15 bold').grid(row=1, column=1,padx=10)
            del print_data[-1]
            for i in range(len(print_data)):
                tkinter.Label(Result, text=round(print_data[i][0], 2), font='Times 15 bold').grid(row=i + 2, column=0,padx=10,pady=5)
                tkinter.Label(Result, text=round(print_data[i][1], 2), font='Times 15 bold').grid(row=i + 2, column=1,padx=10,pady=5)
                tkinter.Label(Result, text=round(print_data[i][2], 2), font='Times 15 bold').grid(row=i + 2, column=2,padx=10,pady=5)
                tkinter.Label(Result, text=print_data[i][4], font='Times 15 bold').grid(row=i + 2, column=3,padx=10,pady=5)
                tkinter.Label(Result, text=print_data[i][3], font='Times 15 bold').grid(row=i + 2, column=4,padx=10,pady=5)
        else:
            tkinter.messagebox.showinfo('Message', '''The  columns named 'Age(Ma)' and 'Error(Ma)' must be contained in the original file''', parent = t4)
    t4 = Toplevel()
    t4.title('Multi-population')
    t4.geometry(s_center)
    t4.geometry('970x805+370+240')
    t4.resizable(True, False)
    Operation = tkinter.LabelFrame(t4, text='Operation', width=250, height=800,font='Times 17 bold')
    Operation.grid(row=0, column=0, padx=5)
    Operation.grid_propagate(False)
    Result = tkinter.LabelFrame(t4, text='Result', width=700, height=800,font='Times 17 bold')
    Result.grid(row=0, column=1)
    Result.grid_propagate(False)
    Data_input = tkinter.LabelFrame(Operation, text='Input', width=240, height=80,font='Times 15 bold')
    Data_input.grid(row=0, column=0, padx=5, pady=5, rowspan=6, columnspan=2, sticky=N)
    Data_input.grid_propagate(False)
    show_file = tkinter.Entry(Data_input, relief='groove', width=18)
    show_file.grid(row=0, column=0, pady=5, padx=5, sticky=W + N + S)
    tkinter.Button(Data_input, text='Open', width=7,font='Times 15', command=secection).grid(row=0, column=5, padx=5, pady=5, columnspan=3,sticky=E)
    Type_error = tkinter.LabelFrame(Operation, text='Type(error)', width=240, height=90,font='Times 15 bold')
    Type_error.grid(row=8, column=0, padx=5, pady=5, rowspan=6, columnspan=2, sticky=N)
    Type_error.grid_propagate(False)
    type_sigma = IntVar()
    type_sigma.set('1')  # The default choice is 2sigma
    Radiobutton(Type_error, variable=type_sigma, value='1', text='1sigma').grid(row=0, column=0, padx=30, pady=10)
    Radiobutton(Type_error, variable=type_sigma, value='2', text='2sigma').grid(row=1, column=0, pady=5)
    type_content = IntVar()
    type_content.set('1')  # The default selection error type is percentage.
    Radiobutton(Type_error, variable=type_content, value='1', text='abs.(ppm)').grid(row=0, column=1, sticky=W)
    Radiobutton(Type_error, variable=type_content, value='2', text='percent(%)').grid(row=1, column=1, sticky=W)
    tkinter.Label(Operation, text='Components:', font='Times 16 bold').grid(row=16, column=0, padx=10, columnspan=2,sticky=W+N+S)
    spin = Spinbox(Operation, from_=1, to=10, width=10)
    spin.insert(0, 1)  # Set the number of initial components to 1.
    spin.grid(row=16, column=0, pady=5, padx=5, columnspan=2, sticky=E)  # sp.get ( ) is used to obtain the current value.
    tkinter.Button(Operation, text='Calculate Age',font='Times 15 ', command=count).grid(row=17, column=0, pady=5, padx=5, columnspan=4,sticky=W + E)
    tkinter.Button(Operation, text='Reset',font='Times 15 ', command=clear_all).grid(row=18, column=0, pady=5, padx=5, columnspan=4,sticky=W + E)
    SSE = tkinter.LabelFrame(Operation, text='Show', width=240, height=93,font='Times 15 bold')
    SSE.grid(row=19, column=0, padx=5, columnspan=2)
    SSE.grid_propagate(False)
    tkinter.Label(SSE, text='P_total:', font='Times 18 bold').grid(row=1, column=0, padx=10, rowspan=1,sticky=S + E)
    tkinter.Label(Result, text='------------------------Calculated------------------------', font='Times 25 bold').grid(row=0, column=0,padx=10,pady=10,columnspan=5,sticky=W + E)
    tkinter.Label(Result, text='Age(Ma)', font='Times 15 bold').grid(row=1, column=0, padx=10, pady=5)
    tkinter.Label(Result, text='Error(Ma)', font='Times 15 bold').grid(row=1, column=1, padx=10, pady=5)
    tkinter.Label(Result, text='MSWD', font='Times 15 bold').grid(row=1, column=2, padx=10, pady=5)
    tkinter.Label(Result, text='P', font='Times 15 bold').grid(row=1, column=3, padx=10, pady=5)
    tkinter.Label(Result, text='N', font='Times 15 bold').grid(row=1, column=4, padx=10, pady=5)

win = tkinter.Tk()
plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams['axes.unicode_minus'] = False
stlyle = ttk.Style(theme='cosmo')
win.title('Copyright © 2023 CAUM. All Rights Reserved')
ctypes.windll.shcore.SetProcessDpiAwareness(1)
width, height = 1400, 1010
width_max, height_max = win.maxsize()
s_center = '%dx%d+%d+%d' % (width, height, (width_max - width) / 0.5, (height_max - height) / 80)
win.geometry(s_center)
win.resizable(False, False)
box1 = tkinter.LabelFrame(win, text='Input', font='Times 15 bold')
box1.grid(row=0, column=0, padx=5, sticky=W + E)
tkinter.Label(box1, text='File path:', height=1, font='Times 16 bold').grid(row=0, column=0, padx=5, pady=10, sticky=S)
path = tkinter.Text(box1, width=61, height=1, relief='groove')
path.grid(row=0, column=1, padx=5, pady=5)
tkinter.Button(box1, text='Open', width=15, font=('Times', 16), command=Open).grid(row=0, column=2, padx=10, pady=10,sticky=N + S, rowspan=2)
tkinter.Button(box1, text='Run', width=15, font=('Times', 16), command=Run).grid(row=0, column=3, padx=10, pady=10,sticky=S + N, rowspan=2)
tkinter.Button(box1, text='Reset', width=15, font=('Times', 16), command=Reset).grid(row=0, column=4, padx=10, pady=10,sticky=S + N, rowspan=2)
tkinter.Button(box1, text='Save Excel', width=15, font=('Times', 16), command=Save_excel).grid(row=0, column=5, padx=10,pady=10, sticky=S + N,rowspan=2)
box_error = tkinter.LabelFrame(box1, text='Error(%)', font='Times 16 bold')
box_error.grid(row=1, column=0, padx=5, pady=15, columnspan=2, sticky=E + W)
tkinter.Label(box_error, text='UO2:', height=1, font='Times 14 bold').grid(row=0, column=0, padx=5, pady=5)
error_UO2 = tkinter.Entry(box_error, relief='groove', width=12)
error_UO2.grid(row=0, column=1, pady=5, sticky=E)
tkinter.Label(box_error, text='%', height=1, font='Times 14 bold').grid(row=0, column=2, padx=0, pady=5)
tkinter.Label(box_error, text='ThO2:', height=1, font='Times 14 bold').grid(row=0, column=3, padx=5, pady=5)
error_ThO2 = tkinter.Entry(box_error, relief='groove', width=12)
error_ThO2.grid(row=0, column=4, pady=5, sticky=E)
tkinter.Label(box_error, text='%', height=1, font='Times 14 bold').grid(row=0, column=5, padx=0, pady=5)
tkinter.Label(box_error, text='PbO:', height=1, font='Times 14 bold').grid(row=0, column=6, padx=5, pady=5)
error_PbO = tkinter.Entry(box_error, relief='groove', width=12)
error_PbO.grid(row=0, column=7, pady=5, sticky=E)
tkinter.Label(box_error, text='%', height=1, font='Times 14 bold').grid(row=0, column=8, padx=0, pady=5)
Age_excel = tkinter.LabelFrame(win, text='-', height=605, font='Times 5 bold', labelanchor='ne')
Age_excel.grid(row=1, column=0, padx=5, pady=10, sticky=W + E)
notebook1 = ttk.Notebook(Age_excel)
frame1 = tkinter.Frame(notebook1, width=615, height=545)
frame2 = tkinter.Frame(notebook1, width=615, height=545)
notebook1.add(frame1, text="The chemical age")
notebook1.add(frame2, text="Coefficient of correlation")
notebook1.grid(row=0, column=0)
frame3 = LabelFrame(Age_excel, text='The probability density plot', width=617, height=561)
frame3.grid(row=0, column=1, padx=5, sticky=S)
# Set the header font
style = ttk.Style()
style.configure("Treeview.Heading", font='Times 15 ')
style.configure('Treeview', font='Times 15 ', rowheight=30)
colnames1 = ['ID', 'UO2(wt%)', 'ThO2(wt%)', 'PbO(wt%)', 'Age(Ma)', 'Error(Ma)']
colnames2 = ['Variable 1', 'Variable 2', 'P', 'R']
tree1 = ttk.Treeview(frame1, columns=colnames1, show='headings', height=17)
tree1.grid(row=0, column=0, sticky=W + E, padx=3)
tree2 = ttk.Treeview(frame2, columns=colnames2, show='headings', height=17)
tree2.grid(row=0, column=0, sticky=W + E + N, padx=3)
# Define the header
tree1.heading('ID', text='ID')
tree1.heading('UO2(wt%)', text='UO2(wt%)')
tree1.heading('ThO2(wt%)', text='ThO2(wt%)')
tree1.heading('PbO(wt%)', text='PbO(wt%)')
tree1.heading('Age(Ma)', text='Age(Ma)')
tree1.heading('Error(Ma)', text='Error(Ma)')
tree2.heading('Variable 1', text='Variable 1')
tree2.heading('Variable 2', text='Variable 2')
tree2.heading('P', text='P')
tree2.heading('R', text='R')
tree1.column('ID', anchor=CENTER, width=80)
tree1.column('UO2(wt%)', anchor=CENTER, width=130)
tree1.column('ThO2(wt%)', anchor=CENTER, width=130)
tree1.column('PbO(wt%)', anchor=CENTER, width=130)
tree1.column('Age(Ma)', anchor=CENTER, width=140)
tree1.column('Error(Ma)', anchor=CENTER, width=140)
tree2.column('Variable 1', anchor=CENTER, width=170)
tree2.column('Variable 2', anchor=CENTER, width=170)
tree2.column('P', anchor=CENTER, width=200)
tree2.column('R', anchor=CENTER, width=210)
# Add a scroll bar on the Y axis
yscrollbar1 = Scrollbar(frame1)
yscrollbar1.grid(row=0, column=0, sticky=N + S + E)
yscrollbar1.config(command=tree1.yview)
tree1.configure(yscrollcommand=yscrollbar1.set)
yscrollbar2 = Scrollbar(frame2)
yscrollbar2.grid(row=0, column=0, sticky=N + S + E)
yscrollbar2.config(command=tree2.yview)
tree2.configure(yscrollcommand=yscrollbar2.set)
box2 = tkinter.LabelFrame(win, text='EMPA Dating', font='Times 15 bold')
box2.grid(row=3, column=0, padx=5, sticky=W + E)
Label(box2, text='Based on the results, please choose one of the following modules to treat the data.',font='Times 27 bold').grid(row=0, column=0, padx=8, pady=10, columnspan=3)
tkinter.Button(box2, text='Single population', width=35, height=4, font=('Times', 16), command=weighted_age).grid(row=1,column=0,padx=15,pady=10)
tkinter.Button(box2, text='Multi-population', width=35, height=4, font=('Times', 16), command=Mix_age).grid(row=1,column=1,pady=10,padx=15)
tkinter.Button(box2, text='Age extrapolation', width=35, height=4, font=('Times', 16), command=estimate_age).grid(row=1,column=2,padx=15,pady=10)
win.mainloop()