import numpy as np

def write_record_in_file2d(text_file,date,bottom_levels,varname,n_bound,n_lev,n_vars,record):
	header = '0 2 957839 %d %d %d 1\n' # n_bound, n_lev, n_vars
	date_string_file = '%Y%m%d %H%M%S\n'
	string0 = '%s -999.0 ' % n_lev
	text_file.write(header % (n_bound,n_lev,n_vars))
	text_file.write(date.strftime(date_string_file))
	text_file.write(' '.join('%d ' % ii for ii in bottom_levels) + '\n')
	text_file.write(varname)
	for bb in range(n_bound):
		text_file.write(string0 + ' '.join(' %.8f' % ii for ii in record[:,bb].tolist()) + '\n')
	return

def write_record_in_file2d_vel(text_file,date,bottom_levels,varname1,varname2,n_bound,n_lev,n_vars,record1,record2):
	header = '0 2 957839 %d %d %d 1\n' # n_bound, n_lev, n_vars
	date_string_file = '%Y%m%d %H%M%S\n'
	string0 = '%s -999.0 ' % n_lev
	text_file.write(header % (n_bound,n_lev,n_vars))
	text_file.write(date.strftime(date_string_file))
	text_file.write(' '.join('%d ' % ii for ii in bottom_levels) + '\n')
	text_file.write(varname1)
	for bb in range(n_bound):
		text_file.write(string0 + ' '.join(' %.8f' % ii for ii in record1[:,bb].tolist()) + '\n')
	text_file.write(varname2)
	for bb in range(n_bound):
		text_file.write(string0 + ' '.join(' %.8f' % ii for ii in record2[:,bb].tolist()) + '\n')
	return


def write_record_in_file1d(text_file,date,varname,n_bound,n_lev,n_vars,record):
	header = '0 2 957839 %d %d %d 1\n' # n_bound, n_lev, n_vars
	date_string_file = '%Y%m%d %H%M%S\n'
	text_file.write(header % (n_bound,n_lev,n_vars))
	text_file.write(date.strftime(date_string_file))
	text_file.write(varname)
	for bb in range(n_bound):
		text_file.write('%.8f \n'% record[bb])
	return
