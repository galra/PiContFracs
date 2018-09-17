import configparser
import json
import re
import os
from postprocfuncs import POSTPROC_FUNCS
from basic_enum_params import BasicEnumPolyParams, IndexedParameterEnumPolyParams


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
                raise
                # exit()


        return self._cast_unknown_param_val(param_val)

    # generates a config filename with a new version
    def _gen_hashtable_filename(self):
        hashtable_file_operation = self.config['hashtable_file_operation']
        filename = self.config['hashtable_file']
        if not filename:
            return

        if not re.match('(.*?)(_v[0-9]+)?(\\.pkl)', filename[:-len('.pkl')]):
            filename = '%s_v1.pkl' % filename[:-len('.pkl')]

        stripped_filename = re.split('_v[0-9]+.pkl$', filename)[0]
        old_versions = [ fn for fn in os.listdir() if re.match('(%s)(_v[0-9.]+)?(\\.pkl)' % stripped_filename, fn) ]
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
            filename = '%s_v%d.pkl' % (filename_base, v)
        else:
            filename = old_versions[-1]
        self.config['hashtable_file'] = filename

    def _gen_postproc_funcs_filter(self):
        self.config['postproc_funcs_filter'] = [ POSTPROC_FUNCS.index(f) for f in self.config['postproc_funcs_filter'] ]

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
        self._gen_postproc_funcs_filter()

        if self.config['hashtable_file_operation'] not in ['generate', 'expand', 'use']:
            raise ValueError('Illegal value for hashtable_file_operation')

    def get_config(self):
        return self.config

    def set(self, section, option, value=None):
        if option == 'postproc_funcs_filter' and value:
            funcs_list = json.loads(value)
            funcs_names_list = [ '"%s"' % POSTPROC_FUNCS[i] for i in funcs_list ]

            value = '[%s]' % (', '.join(funcs_names_list))
        return super().set(section=section, option=option, value=value)

    @staticmethod
    def poly_type_parser(s):
        s = ConfigParser.string_parameter_parser(s)
        if s not in AB_POLYS_TYPES:
            raise ValueError('%s is not a valid type of a polynom' % s)
        return AB_POLYS_TYPES[s]


CONFIG_PARAMS_TYPES = {'poly_coeffs_range':  json.loads,
                       'ulcd_range': json.loads,
                       'const': ConfigParser.string_parameter_parser,
                       'ab_polys_type': ConfigParser.poly_type_parser,
                       'a_poly_size': int,
                       'b_poly_size': int,
                       'a_interlace': int,
                       'b_interlace': int,
                       'print_surprising_nonexp_contfracs': ConfigParser.string_parameter_parser,
                       'gen_hashtable_only': ConfigParser.string_parameter_parser,
                       # json is used to support interlace lists
                       'a_coeffs_range': json.loads,
                       'b_coeffs_range': json.loads,
                       'lhs_type': ConfigParser.string_parameter_parser,
                       'lhs_params': json.loads,
                       'postproc_funcs_filter': json.loads,
                       'i': int,
                       'hashtable_file_operation': ConfigParser.string_parameter_parser,
                       'hashtable_file': ConfigParser.string_parameter_parser,
                       'hashtable_num_of_iterations': int,
                       'hashtable_precision': int}

AB_POLYS_TYPES = {str(BasicEnumPolyParams): BasicEnumPolyParams,
                  str(IndexedParameterEnumPolyParams): IndexedParameterEnumPolyParams}
