import sys
import dill as pickle
import os
import progressbar


from enum_poly_params import BasicEnumPolyParams

does_have_int_roots = BasicEnumPolyParams._does_have_integer_roots
mitm = None

def main(hashtable_filepath):
    global mitm
    with open(hashtable_filepath, 'rb') as input_file:
        mitm = pickle.load(input_file)
    dec_hashtable_keys = list(mitm.dec_hashtable.keys())
    dec_hashtable_keys.remove('parameters')

    for k in progressbar.progressbar(dec_hashtable_keys):
        removed_params_counter = 0
        for i, params in enumerate(mitm.dec_hashtable[k][:]):
            b_polys = params[0][1]
            for b_poly in b_polys:
                if does_have_int_roots(b_poly):
                    mitm.dec_hashtable[k].pop(i - removed_params_counter)
                    removed_params_counter += 1
                    break
        if len(mitm.dec_hashtable[k]) == 0:
            mitm.dec_hashtable.pop(k)

    output_filepath = list(os.path.splitext(hashtable_filepath))
    output_filepath[0] += '_fixed'
    output_filepath = ''.join(output_filepath)
    with open(output_filepath, 'wb') as output_file:
        pickle.dump(mitm, output_file, protocol=pickle.HIGHEST_PROTOCOL)


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('Usage: %s <hashtable_file_path>' % sys.argv[0])
        exit()
    main(sys.argv[1])