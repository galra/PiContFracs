from decimal import Decimal as dec

class DecimalHashTable(dict):
    def __init__(self, accuracy):
        # +1 for the decimal point
        self.accuracy = accuracy + 1

    def _manipulate_key(self, key):
        if not isinstance(key, dec) and not isinstance(key, str):
            raise TypeError('Only Decimal is supported')
        if isinstance(key, str):
            return key
        return key.to_eng_string()[:self.accuracy+1]

    def __setitem__(self, key, value):
        key = self._manipulate_key(key)
        return super().__setitem__(key, value)

    def __getitem__(self, item):
        item = self._manipulate_key(item)
        return super().__getitem__(item)

    def __delitem__(self, key):
        key = self._manipulate_key(key)
        return super().__getitem__(key)

    def __contains__(self, item):
        item = self._manipulate_key(item)
        return super().__contains__(item)