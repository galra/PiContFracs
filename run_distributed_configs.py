import os
import sys
import time

def main(configs_dir):
    files_list = [ os.path.join(configs_dir, f) for f in os.listdir(configs_dir) if os.path.splitext(f)[1] == '.ini' ]
    for f in files_list:
        os.system('start cmd.exe /k python main.py %s' % f)

if __name__ == '__main__':
    if len(sys.argv) == 2:
        main(sys.argv[1])
    else:
        print('%s <configs_dir>' % sys.argv[0])