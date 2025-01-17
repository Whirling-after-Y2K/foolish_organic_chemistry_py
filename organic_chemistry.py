CHEMISTRY_BOND_DICT = {'C': 4, 'O': 2, 'N': 3}
FEATURE_TABLE = {'C': 0b1000000,'O':0b0100000}


class Element:
    def __init__(self, name):
        if type(name) is str:
            name = name.upper()
        else:
            raise TypeError
        self.name = name
        self.bond_list = [None for _ in range(CHEMISTRY_BOND_DICT[name])]
        self.feature = 0b00000
        self.feature += FEATURE_TABLE[name]

    # def __eq__(self, other):
    #     return self.feature == other.feature

    def add_bond(self, target_element, connect_num=1, is_first=True):
        for i in range(connect_num):
            if None in self.bond_list:
                available_pos = self.bond_list.index(None)
            else:
                raise OverflowError
            self.bond_list[available_pos] = target_element
        if is_first:
            target_element.add_bond(self, connect_num, False)
        # self.feature += FEATURE_TABLE[str(connect_num) + target_element]


c1 = Element("C")
c2 = Element("o")
c1.add_bond(c2, 2)

print(c1.bond_list, c1.name)
print(c2.bond_list, c2.name)
