import json

CHEMISTRY_BOND_DICT = {'C': 4, 'O': 2, 'N': 3}
HASH_FEATURE_TABLE = {}

with open("chemistry_feature.json", 'r', encoding='utf-8') as f:
    raw_feature_table = json.loads(f.read())
del f

MAX_FEATURE_NUM = 5  # 前五位用于储存self_name
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
        self.feature = add_feature(self.feature, "self_name", self.name)
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

            # 待实现: 实际调用
            # self.feature = update_feature(self.feature, "direct_connection", target_element, connect_num)

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


def inquire_feature(element_feature: int, feature_class: str, feature: str):
    if feature_class == "self_name":
        return element_feature & (1 << 5)
    return element_feature & (1 << HASH_FEATURE_TABLE[feature_class][feature])


def add_feature(element_feature: int, feature_class: str, feature: str):
    if feature_class == "self_name":
        element_feature |= HASH_FEATURE_TABLE[feature_class][feature]
    element_feature |= (1 << HASH_FEATURE_TABLE[feature_class][feature])
    return element_feature


def del_feature(element_feature: int, feature_class: str, feature: str):
    if feature_class == "self_name":
        raise TypeError
    element_feature &= ~(1 << HASH_FEATURE_TABLE[feature_class][feature])
    return element_feature


def update_feature(element_feature: int, feature_class: str, target_element: Element, connect_num: int):
    feature_name = HASH_FEATURE_TABLE["self_name"][element_feature & (1 << 5)] + "-" + str(
        connect_num) + target_element.name

    return


def connect(target_element_list, is_cyclization=False):
    for i in range(len(target_element_list) - 1):
        target_element_list[i].add_bond(target_element_list[i + 1])
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

for i in c1.bond_list:
    if type(i) is Element:
        print(i.name, i.feature)
print()
for i in c2.bond_list:
    if type(i) is Element:
        print(i.name, i.feature)
    elif type(i) is list:
        for j in i:
            print(j.name, j.feature)
print(c1.bond_list, c1)
print(c2.bond_list, c2.name)
