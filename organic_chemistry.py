CHEMISTRY_BOND_DICT = {'C': 4, 'O': 2, 'N': 3}
FEATURE_TABLE = {'C': 0b1000000, 'O': 0b0100000}


def connect(target_element_list, is_cyclization=False):
    for i in range(len(target_element_list) - 1):
        target_element_list[i].add_bond(target_element_list[i + 1])
    if is_cyclization:
        target_element_list[0].add_bond(target_element_list[-1])


class Element:
    def __init__(self, name):
        if type(name) is str:
            name = name.upper()
        else:
            raise TypeError
        self.name = name
        self.bond_list = [None for _ in range(CHEMISTRY_BOND_DICT[name])]
        self.left_index = 0
        # self.feature = 0b00000
        # self.feature += FEATURE_TABLE[name]

    # def __eq__(self, other):
    #     return self.feature == other.feature

    def add_bond(self, target_element, connect_num=1, is_first=True):
        for i in range(connect_num):
            if self.left_index < CHEMISTRY_BOND_DICT[self.name]:
                self.bond_list[self.left_index] = target_element
                self.left_index += 1
            else:
                raise OverflowError
        if is_first:
            target_element.add_bond(self, connect_num, False)
        # self.feature += FEATURE_TABLE[str(connect_num) + target_element]

    def add_pi_pond(self, target_element_list, is_first=True):  # 注意target_element_list包含自身
        tmp_list = target_element_list.copy()
        if self.left_index < CHEMISTRY_BOND_DICT[self.name]:
            self.bond_list[self.left_index] = tmp_list
            self.left_index += 1
        else:
            raise OverflowError
        if is_first:
            target_element_list.remove(self)
            for target_element in target_element_list:
                target_element.add_pi_pond(tmp_list, False)


c1 = Element("C")
c2 = Element("c")
c3 = Element("c")
c4 = Element("c")
c5 = Element("c")
c6 = Element("c")
c7 = Element("c")
o1 = Element("o")
o2 = Element("o")
o3 = Element("o")
c1.add_bond(o1, 2)
c1.add_bond(c2)
connect([c2, c3, c4, c5, c6, c7], True)

c2.add_pi_pond([c2, c3, c4, c5, c6, c7])

print(c1.bond_list, c1)
print(c2.bond_list, c2.name)
