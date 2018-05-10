import configparser
import json
import re
import os


class ConfigParser(configparser.ConfigParser):
    def __init__(self, configfile=None, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.config = {}
        if configfile:
            self._load_config(configfile)

    @staticmethod
    def string_parameter_parser(s):
        """returns a stripped string, converts None, True and False to the appropriate object"""
        s = s.strip()
        if s.lower() == 'none':
            return None
        elif s.lower() == 'false':
            return False
        elif s.lower() == 'true':
            return True

        return s

    # Cast a string to its optimal cast (int, float, list or still a string)
    def _cast_unknown_param_val(self, s):
        s = self.string_parameter_parser(s)
        if not isinstance(s, str):
            return s

        # Funnel down through integers, floats, lists of numbers and later strings
        try:
            return int(s)
        except ValueError:
            pass

        try:
            return float(s)
        except ValueError:
            pass

        try:
            return json.loads(s)
        except json.JSONDecodeError:
            pass

        return s

    def _cast_param_val(self, param_key, param_val):
        if param_key in CONFIG_PARAMS_TYPES:
            try:
                if param_val.lower().strip() == 'none':
                    return None
                return CONFIG_PARAMS_TYPES[param_key](param_val)
            except json.JSONDecodeError:
                print('Error in _cast_param_val, configfile.py')
                print(param_key)
                print(param_val)
                exit()


        return self._cast_unknown_param_val(param_val)

    # generates a config filename with a new version
    def _gen_hashtable_filename(self):
        hashtable_file_operation = self.config['hashtable_file_operation']
        filename = self.config['hashtable_file']
        if not filename:
            return

        if not re.match('(.*?)(_v[0-9.]+)?(\\.pkl)', filename):
            filename = '%s_v1.pkl' % filename[:-len('.pkl')]

        old_versions = [ fn for fn in os.listdir() if re.match('(.*?)(_v[0-9.]+)?(\\.pkl)', fn) ]
        old_versions.sort()
        if not old_versions:
            self.config['hashtable_file'] = filename
            return
        if hashtable_file_operation == 'generate':
            filename_base, v, _ = re.findall('(.*?)(_v[0-9.]+)?(\\.pkl)', old_versions[-1])[0]
            v = v.replace('_v', '')
            if not v:
                v = 0
            v = int(v) + 1
            print(v)
            filename = '%s_v%d.pkl' % (filename_base, v)
        else:
            filename = old_versions[-1]
        self.config['hashtable_file'] = filename

    # Load a configuration file
    def _load_config(self, file_name='config.ini'):
        self.config = {}

        parsed_config = configparser.ConfigParser()
        results = parsed_config.read(file_name)

        # If no config file was load
        if len(results) == 0:
            return None

        for section in parsed_config.sections():
            for item_key, item_val in parsed_config[section].items():
                self.config[item_key] = self._cast_param_val(item_key, item_val)
        self._gen_hashtable_filename()

        if self.config['hashtable_file_operation'] not in ['generate', 'expand', 'use']:
            raise ValueError('Illegal value for hashtable_file_operation')

    def get_config(self):
        return self.config



CONFIG_PARAMS_TYPES = {'poly_coeffs_range':  json.loads,
                       'ulcd_range': json.loads,
                       'const': ConfigParser.string_parameter_parser,
                       'a_poly_size': int,
                       'b_poly_size': int,
                       'a_interlace': int,
                       'b_interlace': int,
                       'print_surprising_nonexp_contfracs': ConfigParser.string_parameter_parser,
                       # json is used to support interlace lists
                       'a_coeffs_range': json.loads,
                       'b_coeffs_range': json.loads,
                       'u_range': json.loads,
                       'l_range': json.loads,
                       'c_range': json.loads,
                       'd_range': json.loads,
                       'i': int,
                       'hashtable_file_operation': ConfigParser.string_parameter_parser,
                       'hashtable_file': ConfigParser.string_parameter_parser}
