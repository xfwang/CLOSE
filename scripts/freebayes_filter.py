import sys

def get_info_map(info):
    info_map = {}
    fields = info.split(';')
    for field in fields:
        if '=' in field:
            idx = field.index('=')
            name = field[:idx]
            value = field[idx+1:]
            info_map[name] = value
        
    return info_map

min_mq = float(sys.argv[1])
min_mqm = float(sys.argv[2])
max_fs = float(sys.argv[3])
max_read_pos = float(sys.argv[4])
min_dp = float(sys.argv[5])

pass_count = 0
filter_count = 0

for line in sys.stdin:
    line = line.rstrip()
    if line.startswith('#'):
        print line
    else:
        fields = line.split('\t')
        if len(fields) > 8:
            filter = fields[6]
            info = fields[7]
            info_map = get_info_map(info)
            type = info_map['TYPE']
            
            # Absence of FS or ReadPosRankSum indicates indel or no reference support.
            # We do not filter in these cases
            # Absence of other values indicates a problem with the call and these are filtered
            
            mq = 0
            if 'MQ' in info_map:
                mq = float(info_map['MQ'])
            
            mqm = 0
            if 'MQM' in info_map:
                mqm = float(info_map['MQM'].split(',')[0])
            
            fs = 0
            if 'FS' in info_map:
                fs = float(info_map['FS'])
                
            read_pos = 0
            if 'ReadPosRankSum' in info_map:
                read_pos = abs(float(info_map['ReadPosRankSum']))
                
            dp = 0
            if 'DP' in info_map:
                dp = float(info_map['DP']) 
            
            if mq > min_mq and mqm > min_mqm and fs < max_fs and read_pos < max_read_pos and dp > min_dp:
                fields[6] = 'PASS'
                pass_count += 1
            else:
                fields[6] = 'FILTER'
                filter_count += 1

            output = '\t'.join(fields)
            
            print output

            
sys.stderr.write('Passed: ' + str(pass_count) + '\n')
sys.stderr.write('Filtered: ' + str(filter_count) + '\n')
sys.stderr.write('Done.\n')

