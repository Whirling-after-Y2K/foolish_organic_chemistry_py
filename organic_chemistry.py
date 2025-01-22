import json

CHEMISTRY_BOND_DICT = {'c': 4, 'o': 2, 'n': 3}
HASH_FEATURE_TABLE = {}

with open("chemistry_feature.json", 'r', encoding='utf-8') as f:
    raw_feature_table = json.loads(f.read())
del f
SAVE_NAME_NUM = 6
MAX_FEATURE_NUM = SAVE_NAME_NUM  # 前6位用于储存self_name
tmp_index = 1

for i in raw_feature_table:
    if i == "comment":
        continue
    elif i == "self_name":
        HASH_FEATURE_TABLE[i] = {}
        tmp_index2 = 0
        for j in raw_feature_table[i]:
            HASH_FEATURE_TABLE[i][j[0]] = tmp_index2
            HASH_FEATURE_TABLE[i][tmp_index2] = j[0]
            tmp_index2 += 1
        del tmp_index2
    else:
        HASH_FEATURE_TABLE[i] = {}
        for j in raw_feature_table[i]:
            MAX_FEATURE_NUM += 1
            HASH_FEATURE_TABLE[i][j[0]] = tmp_index
            tmp_index += 1

del tmp_index
del raw_feature_table
print(MAX_FEATURE_NUM, HASH_FEATURE_TABLE)


class Element:
    def __init__(self, name: str):
        name = name.lower()
        self.name = name
        self.bond_list = [None for _ in range(CHEMISTRY_BOND_DICT[name])]
        self.feature = 1 << MAX_FEATURE_NUM
        add_feature(self, "self_name", self.name)
        # self.feature = 0b00000
        # self.feature += FEATURE_TABLE[name]

    def add_bond(self, target_element, connect_num=1, is_first=True):
        try:
            left_index = self.bond_list.index(None)
        except ValueError:
            raise IndexError("list index out of range")
        else:
            for _ in range(connect_num):
                self.bond_list[left_index] = target_element
                left_index += 1

            update_feature(self, target_element, connect_num)

            if is_first:
                target_element.add_bond(self, connect_num, False)
        # self.feature += FEATURE_TABLE[str(connect_num) + target_element]

    # def __eq__(self, other):
    #     return self.feature == other.feature

    def add_pi_pond(self, target_element_list: list, is_first=True):  # 注意target_element_list包含自身
        tmp_list: list = target_element_list.copy()
        try:
            left_index: int = self.bond_list.index(None)
        except ValueError:
            raise IndexError("list index out of range")
        else:
            self.bond_list[left_index] = tmp_list
            if is_first:
                target_element_list.remove(self)
                for target_element in target_element_list:
                    target_element.add_pi_pond(tmp_list, False)


def inquire_feature(element: Element, feature_class: str, feature: str):
    element_feature = element.feature
    if feature_class == "self_name":
        return element_feature & (1 << SAVE_NAME_NUM)
    return element_feature & (1 << HASH_FEATURE_TABLE[feature_class][feature])


def add_feature(element: Element, feature_class: str, feature: str):
    if feature_class == "self_name":
        element.feature |= HASH_FEATURE_TABLE[feature_class][feature]
    else:
        element.feature |= (1 << HASH_FEATURE_TABLE[feature_class][feature])


def del_feature(element: Element, feature_class: str, feature: str):
    if feature_class == "self_name":
        raise TypeError
    else:
        element.feature &= ~(1 << HASH_FEATURE_TABLE[feature_class][feature])


def update_feature(element: Element, target_element: Element, connect_num: int, is_first=True):
    if is_first:
        feature_name = element.name + "-" + str(connect_num) + target_element.name
        add_feature(element, "direct_connection", feature_name)
    is_occupy = is_first
    for indirect_connected in target_element.bond_list:
        if indirect_connected is None:
            if is_occupy:
                is_occupy = False
                continue


def connect(target_element_list, is_cyclization=False):
    for tmp in range(len(target_element_list) - 1):
        target_element_list[tmp].add_bond(target_element_list[tmp + 1])
    if is_cyclization:
        target_element_list[0].add_bond(target_element_list[-1])


c1 = Element("C")
c2 = Element("c")
c3 = Element("c")
c4 = Element("c")
c5 = Element("c")
c6 = Element("c")
c7 = Element("c")
# c1-c7 o1 构成苯醛
o1 = Element("o")
o2 = Element("o")
o3 = Element("o")
c1.add_bond(o1, 2)
c1.add_bond(c2)
connect([c2, c3, c4, c5, c6, c7], True)

c2.add_pi_pond([c2, c3, c4, c5, c6, c7])
print(c1.name, id(c1))
print()
for i in c1.bond_list:
    if type(i) is Element:
        print(i.name, i.feature, id(i))
print()
for i in c2.bond_list:
    if type(i) is Element:
        print(i.name, i.feature, id(i))
    elif type(i) is list:
        print()
        for j in i:
            print(j.name, j.feature, id(j))
print(c1.bond_list, c1)
print(c2.bond_list, c2.name)
