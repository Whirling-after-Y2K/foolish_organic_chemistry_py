import json

CHEMISTRY_BOND_DICT = {'c': 4, 'o': 2, 'n': 3}
HASH_FEATURE_TABLE = {}
SAVE_NAME_NUM = 4
MAX_FEATURE_NUM = SAVE_NAME_NUM  # 后4位用于储存self_name 元素名
# MAX_ID_NUM = 11  # 第一位弃用,二到十二位用于储存element_id 最多容纳2^11-2个原子
# element_id = -1

with open("chemistry_feature.json", 'r', encoding='utf-8') as f:
    raw_feature_table = json.loads(f.read())
del f

for i in raw_feature_table:
    if i == "comment":
        continue
    elif i == "self_name":

        tmp_index2 = 0
        for j in raw_feature_table[i]:
            HASH_FEATURE_TABLE[j[0]] = tmp_index2  # tmp_index2用于为每一种名称编号 表示储存在后5位中该元素为tmp_index2
            HASH_FEATURE_TABLE[tmp_index2] = j[0]
            tmp_index2 += 1
        del tmp_index2
    else:
        tmp_index = SAVE_NAME_NUM
        tmp_feature_num = SAVE_NAME_NUM
        for j in raw_feature_table[i]:
            tmp_feature_num += 1  # 每多一种特征就需要多一位进行储存,但不同元素的特征不可能共存 tmp_feature_num用于比较至少需要多少位进行储存
            HASH_FEATURE_TABLE[j[0]] = tmp_index  # tmp_index用于为每一种特征编号 表示其储存在第tmp_index+1位
            tmp_index += 1

        MAX_FEATURE_NUM = max(tmp_feature_num - 1, MAX_FEATURE_NUM)

del tmp_index
del raw_feature_table


# print(MAX_FEATURE_NUM, HASH_FEATURE_TABLE)


class Molecule:
    def __init__(self, name):
        self.name = name
        self.composition = []
        self.feature = []
        # self.feature_dict = {}
        # self.link_dict = {}
        self.c_num = 0
        self.o_num = 0
        self.n_num = 0

    def __eq__(self, other):
        return self.feature == other.feature

    def __hash__(self):
        return hash(self.feature)

    def update_name(self, new_name):
        self.name = new_name

    def update(self):
        self.feature.clear()
        for atom in self.composition:
            atom.overall_f = update_overall_f(atom)
            self.feature.append(atom.overall_f)
        self.feature.sort()


class Atom:
    def __init__(self, name: str, molecule: Molecule):
        name = name.lower()
        self.name = name
        self.bond_list = [None for _ in range(CHEMISTRY_BOND_DICT[name])]
        self.feature = 1 << MAX_FEATURE_NUM
        self.feature |= HASH_FEATURE_TABLE[self.name]
        self.overall_f = ''
        self.belong = molecule
        molecule.composition.append(self)

    def add_bond(self, target_atom, connect_num=1, is_first=True):
        feature_name = self.name + "-" + str(connect_num) + target_atom.name
        add_feature(self, feature_name)
        if is_first:
            target_atom.add_bond(self, connect_num, False)

        add_point = self.bond_list.index(None)
        for _ in range(connect_num):
            self.bond_list[add_point] = target_atom
            add_point += 1

    def add_pi_pond(self, target_atom_list: list, is_first=True):  # 注意target_element_list包含自身
        tmp_list: list = target_atom_list.copy()
        self.bond_list[self.bond_list.index(None)] = tmp_list
        if is_first:
            target_atom_list.remove(self)
            for target_atom in target_atom_list:
                target_atom.add_pi_pond(tmp_list, False)


'''
def inquire_id(atom_feature: int):
    return (atom_feature & (((1 << (MAX_FEATURE_NUM + MAX_ID_NUM)) - 1) & (~((1 << MAX_FEATURE_NUM) - 1))))>>MAX_FEATURE_NUM
# 用于实现对id的查询,至于为什么不用字符串切片,大概是因为位运算它快吧
# 此函数已正式弃用,但因为位运算很帅,所以将其保留
'''


def inquire_feature(atom_feature: int, feature: str):
    return atom_feature & (1 << HASH_FEATURE_TABLE[feature])


def add_feature(atom: Atom, feature: str):
    atom.feature |= (1 << HASH_FEATURE_TABLE[feature])


def del_feature(atom: Atom, feature: str):
    atom.feature &= ~(1 << HASH_FEATURE_TABLE[feature])


# update_overall_f函数,通过BFS实现从自己开始逐步获取所连原子的feature属性
def update_overall_f(self):
    visit_point = 0
    will_visit = [self]
    distant = 0
    # self.overall_f += father_f
    ansL = [str(distant) + str(self.feature)]

    while visit_point < len(will_visit):
        self = will_visit[visit_point]
        distant = int(ansL[visit_point][0]) + 1
        for fa in self.bond_list:
            if (fa is None) or (fa in will_visit) or (type(fa) is list):
                continue
            else:
                will_visit.append(fa)
                ansL.append(str(distant) + str(fa.feature))
        visit_point += 1
    ansL.sort()
    return '.'.join(ansL)


def connect(target_atom_list, is_cyclization=False):
    for tmp in range(len(target_atom_list) - 1):
        target_atom_list[tmp].add_bond(target_atom_list[tmp + 1])
    if is_cyclization:
        target_atom_list[0].add_bond(target_atom_list[-1])


if __name__ == '__main__':
    benzaldehyde = Molecule("benzaldehyde")
    # benzaldehyde.add_atom('c', 7)
    # benzaldehyde.add_atom('o')

    c1 = Atom("C", benzaldehyde)
    c2 = Atom("c", benzaldehyde)
    c3 = Atom("c", benzaldehyde)
    c4 = Atom("c", benzaldehyde)
    c5 = Atom("c", benzaldehyde)
    c6 = Atom("c", benzaldehyde)
    c7 = Atom("c", benzaldehyde)
    # c1-c7 o1 构成苯醛
    o1 = Atom("o", benzaldehyde)
    # o2 = Atom("o", benzaldehyde)
    # o3 = Atom("o", benzaldehyde)

    c1.add_bond(o1, 2)
    c1.add_bond(c2)
    connect([c2, c3, c4, c5, c6, c7], True)

    benzaldehyde1 = Molecule("benzaldehyde1")
    # benzaldehyde1.add_atom('c', 7)
    # benzaldehyde1.add_atom('o')

    c11 = Atom("C", benzaldehyde1)
    c21 = Atom("c", benzaldehyde1)
    c31 = Atom("c", benzaldehyde1)
    c41 = Atom("c", benzaldehyde1)
    c51 = Atom("c", benzaldehyde1)
    c611 = Atom("c", benzaldehyde1)
    c71 = Atom("c", benzaldehyde1)
    # c1-c7 o1 构成苯醛
    o11 = Atom("o", benzaldehyde1)
    # o2 = Atom("o", benzaldehyde1)
    # o3 = Atom("o", benzaldehyde1)

    c11.add_bond(o11, 2)
    c11.add_bond(c21)
    connect([c611, c41, c31, c51, c21, c71], True)

    benzaldehyde.update()
    benzaldehyde1.update()
    print(benzaldehyde.feature)
    print(bool(benzaldehyde == benzaldehyde1))
    for i in c1.bond_list:
        if i is None:
            continue
        print(i.overall_f)

    c2.add_pi_pond([c2, c3, c4, c5, c6, c7])
    print(c1.name, c1.overall_f)
    print()
    print(benzaldehyde.composition)
    for i in benzaldehyde.composition:
        print(i.name, i.feature, i.overall_f)

    # print(c1.bond_list, c1)
    # print(c2.bond_list, c2.name)

