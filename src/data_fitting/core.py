import os
import colorsys
import matplotlib
import numpy as np
import pandas as pd
from inspect import signature
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import matplotlib.font_manager as fm
from matplotlib.ticker import FormatStrFormatter

# Set font embedding type for Illustrator compatibility
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

# Try to load Myriad Pro, fallback to default
MYRIAD_FONT_PATH = "/Library/Fonts/Myriad Pro Regular.ttf"
if os.path.exists(MYRIAD_FONT_PATH):
    default_prop = fm.FontProperties(fname=MYRIAD_FONT_PATH, size=24)
    legend_prop = fm.FontProperties(fname=MYRIAD_FONT_PATH, size=16)
else:
    default_prop = fm.FontProperties(size=24)
    legend_prop = fm.FontProperties(size=16)

def logistic_HS(c, T, B, HS, EC50):
	'''
	logistic equation used primarily for binding curves. 
	This is a 4-parameter equation with vairable Hill Slope
	'''
	return c**HS * (T - B) / (c ** HS + EC50 ** HS) + B

def logistic(c, T, B, EC50):
	'''
	logistic equation used primarily for binding curves. 
	This is a 3-parameter equation with unit Hill Slope
	'''
	return c * (T-B) / (c + EC50) + B

def bi_modal(T, B, p, mu1, s1, Fmax1, mu2, s2, Fmax2):
    phi1 = Fmax1/(s1 * np.sqrt(2 * np.pi)) * np.exp(-1/2 * ((T - mu1) / s1) ** (2))
    phi2 = Fmax2/(s2 * np.sqrt(2 * np.pi)) * np.exp(-1/2 * ((T - mu2) / s2) ** (2))
    
    return B + p * phi1 + (1 - p) * phi2

def get_fitting_parameters(curve_function):
	'''
	get_fitting_parameters gets the names of the parameters of the curve you are trying to fit.
	input:
		curve_function:
			<class 'function'> function of the curve you are trying to fit data to.
			curve_function must take independent variable as the first argument and the remainig arguments 
			are the parameters of the model.
	output:
		params:
			<class 'list'> list of the parameter names of curve_function
		output:
			<class 'dict'> dictionairy of fitted parameters
		stdev:
			<class 'dict'> dictionary of standard deviation in fitted parameter
	'''
	params = str(signature(curve_function)).replace("(",'').replace(")",'').replace(' ', '').split(',')[1::]
	output = {}
	stdev = {}
	for i in range(len(params)):
		output[params[i]] = []
		stdev['{} stdev'.format(params[i])] = []
	
	return params, output, stdev

def get_data(data_excel, sheet_name = 'RAW DATA', error_bars = False):
	'''
	get_data retrieves the data to be fitted from data_excel.

	input:
		data_excel:
			<class 'str'> string of excel file with data to fit. if the file is not in 
			your current directory the path to the file will also work
			the excel file must contain the data in a specified format.
				
				'RAW DATA' tab contains the data to fit. the first column
				is the independent variable and the remaining columns are the
				dependent measurements. These columns must contain the same 
				number of values. if any measurement is missing a data point
				it will need to be in its own excel sheet. the headder of each
				measurement column will be used for the legend on the plot.

				'AXES LABEL' tab contains the x-axis label in cell A1 and 
				the y-axis label in cell A2.

				'ERROR BAR' tab is optional and used to specify the error 
				bars for repeated measurements. This tab is usually used
				for plotting the averaged data and not for obtaining the 
				fitted parameters .

	output:
		xdata:
			<class 'list'> of <class 'float'> list of independent variable values from data_excel:RAW_DATA
		ydata:
			<class 'list'> of <class 'list'> of <class 'float'> list of list of each dependent variable values from data_excel:RAW_DATA
		legends:
			<class 'list'> list of dependent variable labels cells A2-?2 one for each set of data to fit
		axes_label:
			<class 'list'> list of x-axis label and y-axis label from data_excel:AXES_LABEL
		y_error:
			<class 'list'> of <class 'listl> similar to ydata but contains the specified error for each measurement
			from data_excel:ERROR BAR
	'''
	#gets raw data from experiment file in sheet named RAW DATA
	df = pd.read_excel(data_excel, sheet_name = sheet_name, index_col = 0)

	#gets axes label infor from AXES LABEL sheet
	df2 = pd.read_excel(data_excel, sheet_name = 'AXES LABEL')
	axes_label = list(df2.columns)

	#gets error bars if they exist.
	if error_bars:
		df3 = pd.read_excel(data_excel, sheet_name = error_bars, index_col = 0)
		y_error = []
		for i in range(len(df3.columns)):
			y_error.append(list(df3[df3.columns[i]]))
	else:
		y_error = None

	#generates x and y data as well as legends list
	xdata = list(df.index)

	ydata = []
	legends = []

	for i in range(0, len(df.columns)):
		ydata.append(list(df[df.columns[i]]))
		legends.append(df.columns[i])
	
	if y_error is not None:
		return xdata, ydata, legends, axes_label, y_error
	else:
		return xdata, ydata, legends, axes_label,

def get_data_with_missing_values(data_excel, error_bars = False, sheet_name = 'RAW DATA'):
	if error_bars:
		xdata, ydata, legends, axes_label, y_error = get_data(data_excel, sheet_name, error_bars)
	else:
		xdata, ydata, legends, axes_label = get_data(data_excel, sheet_name, error_bars)
		y_error = None
	
	fixed_data = {d:[] for d in legends}

	for j, dataset in enumerate(ydata):
		for i, d in enumerate(dataset):
			if not pd.isna(d):
				try:
					fixed_data[legends[j]][0].append(xdata[i])
					fixed_data[legends[j]][1].append(d)
					if error_bars:
						fixed_data[legends[j]][2].append(y_error[j][i])
				except:
					fixed_data[legends[j]] = [[xdata[i]], [d]]
					if error_bars:
						fixed_data[legends[j]].append([y_error[j][i]])

	return fixed_data, legends, axes_label

def fit_data(
	data_excel, curve_function, number_fit_points = 100, sheet_name = 'RAW DATA', graph = True, bounds = None, 
	log_x = False, log_y = False, fig_size = None, error_bars = False, markersize = 15, save_fig = None,
	colors = None, face = None, xlims = None, ylims = None):

	params, output, stdev = get_fitting_parameters(curve_function)
	print(f'ERROR BARS {error_bars}')
	if graph:
		if fig_size is not None:
			plt.rcParams.update({'font.size': 22})
			fig = plt.figure(figsize = fig_size)
			ax = fig.add_subplot(111)
			
		else:
			fig = plt.figure()
			ax = fig.add_subplot(111)

	if error_bars:
		fixed_data, legends, axes_label = get_data_with_missing_values(data_excel, sheet_name = sheet_name, error_bars = error_bars)

	else:
		print('No Error Bars')
		fixed_data, legends, axes_label = get_data_with_missing_values(data_excel, sheet_name = sheet_name, error_bars = error_bars)
		y_error = None


	for i, (l, data) in enumerate(fixed_data.items()):
		if bounds is not None:
			popt, pcov = curve_fit(curve_function, data[0], data[1], bounds = bounds, maxfev=50000)
		else:
			popt, pcov = curve_fit(curve_function, data[0], data[1], maxfev=50000)

		for j, p in enumerate(popt):
			output[params[j]].append(p)
			stdev[f'{params[j]} stdev'].append(np.sqrt(pcov[j,j])) #might need to verify that this is actually stdev...

		if graph:
			if error_bars:
				y_error = data[2]
			graph_data(ax, data[0], data[1], popt, legends, axes_label, curve_function, number_fit_points, i, log_x, log_y, colors, face, y_error, markersize, xlims, ylims)

	if save_fig is not None:
		plt.savefig(save_fig, dpi = 600, bbox_inches='tight')

	if graph:
		plt.show()

	o_s = {}
	for key in output.keys():
		o_s[key] = []
		for i in range(len(output[key])):
			o_s[key].append('{} +/- {}'.format(round(output[key][i],3),round(stdev['{} stdev'.format(key)][i],3)))

	return output, stdev


def get_bar_graph_data(data_excel, sheet_name = 'BAR GRAPH'):
	'''
	'''
	#gets raw data from experiment file in sheet named RAW DATA
	df = pd.read_excel(data_excel, sheet_name = sheet_name, index_col = 0)

	#gets axes label infor from AXES LABEL sheet
	df2 = pd.read_excel(data_excel, sheet_name = 'AXES LABEL')
	axes_label = list(df2.columns)


def graph_data(axes, xdata, ydata, curve_function_parameters, legends, axes_label, curve_function, num_points, i, log_x, log_y, colors, face, y_error = None, markersize=15, xlims = None, ylims = None):
	'''
	graph_data produces a graph of the data points entered in data_excel:RAW DATA with the fitted curve overlayed in the
	same color.

	input:
		axes:
			<class 'matplotlib.axes._subplots.AxesSubplot'> ax to plot data on
		xdata:
			<class 'list'> xdata list
		ydata:
			<class 'list'> ydata lists one for each set of ydata to plot
		curve_function_parameters:
			<class 'list'> parameter output from scipy.curve_fit
		legends:
			<class 'list'> list of names for each data being fit
		axes_label:
			<class 'list'> list of axis labels
		curve_function:
			<class 'function.> function that data is being fit to
		num_points:
			<class 'int'> number of points to use in the fitted curve. 
			higher numbers will result in smoother curves
		i:
			<class 'int'> current index of which data is being plotted
		log_x:
			<class 'bool'> if True x-axis will be log scale.
		colors:
			<class 'list'> of color identifiers for plotting. 
			for small number of curves being fit, < 11, 
			[plt.get_cmap('tab10')(i) for i in range(10)] works well
		face:
			<class 'list'> of colors for marker face, sometimes no face color
			is wanted so face will look like
			['red', 'none', 'blue', 'none', ...]
		y_error = None:
			<class 'list'> of list of standard deviations for each dependent 
			variable dataset.
		markersize = 15:
			<class 'int'> marker size in the graph.
		xlims = None:
			<class 'list'> [x_low, x_high] used to manually set x limits
		ylims = None:
			<class 'list'> [y_low, y_high] used to manually set y limits
	'''
	import matplotlib
	fpath = '/Library/Fonts/Myriad Pro Regular.ttf'
	if os.path.exists(fpath):
		prop = fm.FontProperties(fname=fpath)
		
	if colors is None:
		colors = [colorsys.hsv_to_rgb(x*1.0/(len(ydata)), 0.75, 0.75) for x in range((len(ydata)))]
		print('set colors')
	if face is None:
		face = [color for color in colors]
		print('set face color')
	xdata_fit = np.linspace(xdata[0], xdata[-1], num_points)
	ydata_fit = []


	for t in xdata_fit:
		ydata_fit.append(curve_function(t, *curve_function_parameters))
	if log_x:
		axes.set_xscale("log")
	if log_y:
		axes.set_yscale('log')
		
	axes.plot(xdata_fit, ydata_fit, color = colors[i], linewidth = 1)
	if y_error is not None:
		axes.errorbar(xdata, ydata,  y_error, fmt = '.', color = colors[i], markersize = markersize, capsize = 2, label = legends[i], markerfacecolor=face[i], ecolor = 'grey')
	else:
		axes.plot(xdata, ydata, '.', color = colors[i], markersize = markersize, label = legends[i], markerfacecolor=face[i])
	
	if xlims is not None:
		axes.set_xlim(xlims)
	if ylims is not None:
		axes.set_ylim(ylims)

	try:
		axes.legend(loc=2, frameon = False, prop = lprop)#, bbox_to_anchor=(1, 0.5),)
	except:
		axes.legend(loc=2, frameon = False)#, bbox_to_anchor=(1, 0.5))

	axes.xaxis.set_major_formatter(FormatStrFormatter('%g'))

	#plt.grid(True)
	try:
		axes.set_xlabel(axes_label[0], fontproperties = prop)
		axes.set_ylabel(axes_label[1], fontproperties = prop)
		for tick in axes.get_xticklabels():
			tick.set_fontproperties(prop)
		for tick in axes.get_yticklabels():
			tick.set_fontproperties(prop)
	except:
		axes.set_xlabel(axes_label[0])
		axes.set_ylabel(axes_label[1])		


