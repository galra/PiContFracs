from decimal import Decimal as dec

# TODO: might be possible to enhance efficiency by using dec.quantize to round to the required
class DecimalHashTable(dict):
    def __init__(self, accuracy):
        # +1 for the decimal point
        self.accuracy = accuracy + 1
        self.accuracy_history = []

    def update_accuracy(self, accuracy):
        self.accuracy_history.append(self.accuracy)
        self.accuracy = accuracy + 1

    def _manipulate_key(self, key):
        if not isinstance(key, dec) and not isinstance(key, str):
            raise TypeError('Only Decimal is supported')
        if isinstance(key, str):
            key_str = key
        else:
            key_str = key.to_eng_string()
        return [ key_str[:i+1] for i in self.accuracy_history ], (key_str[:self.accuracy+1])

    def __setitem__(self, key, value):
        old_keys, cur_key = self._manipulate_key(key)
        for k in old_keys:
            if super().__contains__(k):
                super().__delitem__(k)
        return super().__setitem__(cur_key, value)

    def __getitem__(self, item):
        old_keys, cur_key = self._manipulate_key(item)
        for k in old_keys:
            if super().__contains__(k):
                return super().__getitem__(k)
        return super().__getitem__(cur_key)

    def __delitem__(self, key):
        old_keys, cur_key = self._manipulate_key(key)
        for k in old_keys:
            if super().__contains__(k):
                return super().__delitem__(k)
        return super().__getitem__(cur_key)

    def __contains__(self, item):
        old_keys, cur_key = self._manipulate_key(item)
        for k in old_keys:
            if super().__contains__(k):
                return True
        return super().__contains__(cur_key)

    def __getstate__(self):
        return (self.accuracy, self.accuracy_history, dict(self))

    def __setstate__(self, state):
        self.accuracy, self.accuracy_history, data = state
        self.update(data)

    def __reduce__(self):
        return (DecimalHashTable, (self.accuracy,), self.__getstate__())

    def append_dict(self, appended_dict):
        if self.accuracy != appended_dict.accuracy:
            raise TypeError('Two dictionaries are of non-fitting accuracies')

        self.accuracy_history += [ a for a in appended_dict.accuracy_history if a not in self.accuracy_history ]
        for k in super(DecimalHashTable, appended_dict).keys():
            if k in super().keys():
                try:
                    orig_items = super().__getitem__(k)
                    appended_items = super(DecimalHashTable, appended_dict).__getitem__(k)
                    if not isinstance(orig_items, list):
                        orig_items = [orig_items]
                    if not isinstance(appended_items, list):
                        appended_items = [appended_items]
                    super().__setitem__(k, orig_items + appended_items)
                except TypeError:
                    type_orig = str(type(super().__getitem__(k)))
                    type_appended = str(type(super(DecimalHashTable, appended_dict).__getitem__(k)))
                    print('types are: original dict: %s, appended dict: %s' % (type_orig, type_appended))
                    raise
            else:
                super().__setitem__(k, super(DecimalHashTable, appended_dict).__getitem__(k))