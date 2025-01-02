import cdsapi

c = cdsapi.Client()
root_folder = '/okyanus/users/milicak/dataset/SSHUGVG/global/'

yearind = 2017
for year in range(yearind, yearind+1):
    filename = root_folder + 'global_SSH_UGBG_' + str(year) + '.tar.gz'
    c.retrieve(
        'satellite-sea-level-global',
        {
            'version': 'vDT2021',
            'variable': 'all',
            'format': 'tgz',
            'year': year,
            'month': [
                '01', '02', '03',
                '04', '05', '06',
                '07', '08', '09',
                '10', '11', '12',
            ],
            'day': [
                '01', '02', '03',
                '04', '05', '06',
                '07', '08', '09',
                '10', '11', '12',
                '13', '14', '15',
                '16', '17', '18',
                '19', '20', '21',
                '22', '23', '24',
                '25', '26', '27',
                '28', '29', '30',
                '31',
            ],
        },
        filename)
