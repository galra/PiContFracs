import sys
import os
import enum_params
import dill as pickle

def main(hashtables_dirname, output_file='joined_hashtable.pkl'):
    if os.path.exists(output_file):
        print('Warning: %s exists already. Aborting.' % output_file)
    files_list = [ os.path.join(hashtables_dirname, f) for f in os.listdir(hashtables_dirname)
                   if os.path.splitext(f)[1] == '.pkl' ]
    first_file = files_list.pop(0)
    with open(first_file, 'rb') as f:
        mitm = pickle.load(f)
    for filename in files_list:
        with open(filename, 'rb') as f:
            appended_mitm = pickle.load(f)
        mitm.dec_hashtable.append_dict(appended_mitm.dec_hashtable)
    with open(output_file, 'wb') as out_file:
        pickle.dump(mitm, out_file, pickle.HIGHEST_PROTOCOL)


if __name__ == '__main__':
    if len(sys.argv) == 2:
        main(sys.argv[1])
    elif len(sys.argv) == 3:
        main(sys.argv[1], sys.argv[2])
    else:
        print('%s <hashtables_dir_path> [output_file]\ndefault output_file is joined_hashtable.pkl' % sys.argv[0])