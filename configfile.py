import configparser
import json
import re
import os
from postprocfuncs import POSTPROC_FUNCS
from basic_enum_params import BasicEnumPolyParams, IndexedParameterEnumPolyParams


class ConfigParser(configparser.ConfigParser):
    """Parses a config file in a *.ini format."""
    # To add a value add it to the CONFIG_PARAMS_TYPES dictionary as a key and a type class as a value.
    # int for integers, json.loads for array/lists, ConfigParser.string_parameter_parser for strings
    # (supports none, true, false) or create a new method (e.g. poly_type_parser)
    def __init__(self, configfile=None, *args, **kwargs):
        """Loads and parses a config file.
        Check CONFIG_PARAMS_TYPES for supported and required parameters, others will be ignored.
        See POSTPROC_FUNCS and postprocfuncs.py for supported postprocfuncs and how to add new ones.
        *args and **kwargs are passed to configparser.ConfigParser, see for documentation.
        In the case of an error (a parameter in a wrong format) an exception will be raised."""
        super().__init__(*args, **kwargs)

        # read file, make sure it's not empty
        self._config = {}
        if not configfile:
            return
        parsed_config = configparser.ConfigParser()
        results = parsed_config.read(configfile)
        # If no config file was load
        if len(results) == 0:
            return

        # parse all the parameters
        for section in parsed_config.sections():
            for item_key, item_val in parsed_config[section].items():
                self._config[item_key] = self._cast_param_val(item_key, item_val)
        # update the hashtable name to the new/a new version.
        self._gen_hashtable_filename()
        # convert the postproc functions to their indices
        self._gen_postproc_funcs_filter()

        if self._config['hashtable_file_operation'] not in ['generate', 'expand', 'use']:
            raise ValueError('Illegal value for hashtable_file_operation')

    def _cast_param_val(self, param_key, param_val):
        """Casts parameters to their proper value/type according to CONFIG_PARAMS_TYPES.
        If a parameter isn't defined in CONFIG_PARAMS_TYPES it's ignored."""
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
        return self._cast_unknown_param_val(param_val)

    def _gen_hashtable_filename(self):
        """Generates a config filename with a new version/the newest version (depending on the operation).
        May raise a ValueError exception if the filename doesn't follow the format filename(_v[xx]).pkl"""
        hashtable_file_operation = self._config['hashtable_file_operation']
        filename = self._config['hashtable_file']
        if not filename:
            return

        # validate filename's extension .pkl is specified, or add it
        if not re.match('(.*?)(\\.pkl)', filename):
            filename += '.pkl'

        # if filename is given without a version indication, set a default version v1
        if not re.match('(.*?)(_v[0-9]+)(\\.pkl)', filename):
            filename = '%s_v1.pkl' % filename[:-len('.pkl')]

        # find the oldest version
        stripped_filename = re.split('_v[0-9]+\\.pkl$', filename)[0]
        old_versions = [ fn for fn in os.listdir() if re.match('(%s)(_v[0-9.]+)?(\\.pkl)' % stripped_filename, fn) ]
        old_versions.sort()
        if not old_versions:
            self._config['hashtable_file'] = filename
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
        self._config['hashtable_file'] = filename

    def _gen_postproc_funcs_filter(self):
        """Converts the postproc functions to their indices in POSTPROC_FUNCS"""
        try:
            self._config['postproc_funcs_filter'] = [POSTPROC_FUNCS.index(f) for f in self._config['postproc_funcs_filter']]
        except KeyError as e:
            print('%s is not a valid postproc function' % e.args[0])
            raise

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

    @staticmethod
    def poly_type_parser(s):
        """Converts the poly_type value to the appropriate class."""
        s = ConfigParser.string_parameter_parser(s)
        if s not in AB_POLYS_TYPES:
            raise ValueError('%s is not a valid type of a polynom' % s)
        return AB_POLYS_TYPES[s]

    def _cast_unknown_param_val(self, s):
        """Cast a string to its optimal cast (int, float, list or still a string)."""
        s = self.string_parameter_parser(s)
        if not isinstance(s, str):
            return s

        # Tunnel down through integers, floats, lists of numbers and later strings
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

    def set(self, section, option, value=None):
        """This function is called when generating a config file.
        It wraps the parent class's set to convert the postproc_funcs from indices back to names."""
        if option == 'postproc_funcs_filter' and value:
            funcs_list = json.loads(value)
            funcs_names_list = [ '"%s"' % POSTPROC_FUNCS[i] for i in funcs_list ]

            value = '[%s]' % (', '.join(funcs_names_list))
        return super().set(section=section, option=option, value=value)

    def get_config(self):
        """Returns the parsed config dictionary."""
        return self._config

# defines all the supported config parameters and their types & parser
CONFIG_PARAMS_TYPES = {'const': ConfigParser.string_parameter_parser,
                       'ab_polys_type': ConfigParser.poly_type_parser,
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

# defines the types of supported polynomials
AB_POLYS_TYPES = {str(BasicEnumPolyParams): BasicEnumPolyParams,
                  str(IndexedParameterEnumPolyParams): IndexedParameterEnumPolyParams}
