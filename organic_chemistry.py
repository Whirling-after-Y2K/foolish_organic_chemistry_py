# import json

# -------初始化

HASH_NAME_TABLE = {}
HASH_FEATURE_TABLE = {}

SAVE_NAME_NUM = 3  # 后3位用于储存self_name 元素名 最多8个
# SAVE_FEATURE_NUM = 1  # 用SAVE_FEATURE_NUM位储存特征的数量
# MAX_FEATURE_NUM = SAVE_NAME_NUM

# with open("chemistry_feature.json", 'r', encoding='utf-8') as f:
#     raw_feature_table = json.loads(f.read())
# del f

# del raw_feature_table
CHEMISTRY_BOND_DICT = {"c": 4, "n": 3, "o": 2}
IGNORE_ELEMENT_LIST = ["h", "f", "cl", "br", "i"]

tmp_index1 = 0
for i in CHEMISTRY_BOND_DICT:
    HASH_NAME_TABLE[i] = tmp_index1  # tmp_index1用于为每一种可计算元素编号 表示储存在后SAVE_NAME_NUM位中该元素为tmp_index1
    HASH_NAME_TABLE[tmp_index1] = i
    tmp_index1 += 1

tmp_index1 = SAVE_NAME_NUM
HASH_FEATURE_TABLE['-1h'] = tmp_index1
tmp_index1 += 2  # 特殊用两位储存-h的数量
# 对于一个特征 用1位进行储存 储存在倒数第tmp_index+1位
for i in CHEMISTRY_BOND_DICT:
    feature_name = '-1' + i
    HASH_FEATURE_TABLE[feature_name] = tmp_index1
    tmp_index1 += 1

    feature_name = '-2' + i
    HASH_FEATURE_TABLE[feature_name] = tmp_index1
    tmp_index1 += 1

    if CHEMISTRY_BOND_DICT[i] >= 3:
        feature_name = '-3' + i
        HASH_FEATURE_TABLE[feature_name] = tmp_index1
        tmp_index1 += 1

for i in IGNORE_ELEMENT_LIST:
    if i == 'h':
        continue
    feature_name = '-1' + i
    HASH_FEATURE_TABLE[feature_name] = tmp_index1
    tmp_index1 += 1

HASH_FEATURE_TABLE['c6-'] = tmp_index1
tmp_index1 += 1
HASH_FEATURE_TABLE['no2-'] = tmp_index1
del tmp_index1

# 每一种特征就需要多特定位进行储存
MAX_FEATURE_NUM = len(HASH_FEATURE_TABLE) + (len(HASH_NAME_TABLE) >> 1)


# -------

class Molecule:
    def __init__(self):
        # self.name = name
        self.composition = []
        self.feature = []
        self.unsaturation = 0
        self.formula = {'br': 0, 'cl': 0, 'f': 0, 'h': 0, 'i': 0, 'c': 0, 'n': 0, 'o': 0}
        # self.feature_dict = {}
        # self.link_dict = {}
        # self.c_num = 0
        # self.o_num = 0
        # self.n_num = 0

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
        self.bond_list = ['-1h' for _ in range(CHEMISTRY_BOND_DICT[name])]
        self.feature = 1 << MAX_FEATURE_NUM
        self.feature |= HASH_NAME_TABLE[self.name]
        self.overall_f = ''
        self.belong = molecule
        molecule.composition.append(self)
        molecule.formula[self.name] += 1

    def remove_atom(self, del_target):
        self.bond_list[self.bond_list.index(del_target)] = '-1h'
        self.feature += (1 << SAVE_NAME_NUM)

    def append_atom(self, add_target):
        self.bond_list[self.bond_list.index('-1h')] = add_target
        self.feature -= (1 << SAVE_NAME_NUM)


def add_bond(target_atom0, target_atom1, connect_num=1):
    connected_num = target_atom0.bond_list.count(target_atom1)
    if connected_num > 0:
        target_atom0.belong.unsaturation += connect_num

        feature_name = "-" + str(connected_num) + target_atom1.name
        del_feature(target_atom0, feature_name)
        feature_name = "-" + str(connected_num + connect_num) + target_atom1.name
    else:
        if connect_num > 1:
            target_atom0.belong.unsaturation += (connect_num - 1)
        feature_name = "-" + str(connect_num) + target_atom1.name
    add_feature(target_atom0, feature_name)
    # print(add_point)
    for _ in range(connect_num):
        target_atom0.append_atom(target_atom1)

    target_atom0, target_atom1 = target_atom1, target_atom0

    connected_num = target_atom0.bond_list.count(target_atom1)
    if connected_num > 0:
        feature_name = "-" + str(connected_num) + target_atom1.name
        del_feature(target_atom0, feature_name)
        feature_name = "-" + str(connected_num + connect_num) + target_atom1.name
    else:
        feature_name = "-" + str(connect_num) + target_atom1.name
    add_feature(target_atom0, feature_name)
    # print(add_point)
    for _ in range(connect_num):
        target_atom0.append_atom(target_atom1)
        # add_point += 1


def break_bond(target_atom0, target_atom1, single=False):
    target_atom0.belong.unsaturation -= (target_atom0.bond_list.count(target_atom1) - 1)

    del_feature(target_atom0, '-' + str(target_atom0.bond_list.count(target_atom1)) + target_atom1.name)
    while target_atom1 in target_atom0.bond_list:
        target_atom0.remove_atom(target_atom1)
        # if single:
        #     break

    target_atom0, target_atom1 = target_atom1, target_atom0

    del_feature(target_atom0, '-' + str(target_atom0.bond_list.count(target_atom1)) + target_atom1.name)
    while target_atom1 in target_atom0.bond_list:
        target_atom0.remove_atom(target_atom1)
        # if single:
        #     break


def add_pi_bond(target_atom_list: list):
    if any(i.name == 'c' for i in target_atom_list):
        feature_name = 'c6-'
        target_atom_list[0].belong.unsaturation += 4
    else:
        feature_name = 'no2-'
        target_atom_list[0].belong.unsaturation += 1
    for target_atom in target_atom_list:
        add_feature(target_atom, feature_name)
        target_atom.append_atom(target_atom_list)


def break_pi_bond(target_atom_list: list):
    if any(i.name == 'c' for i in target_atom_list):
        feature_name = 'c6-'
        target_atom_list[0].belong.unsaturation -= 4
    else:
        feature_name = 'no2-'
        target_atom_list[0].belong.unsaturation -= 1
    for target_atom in target_atom_list:
        del_feature(target_atom, feature_name)
        target_atom.remove_atom(target_atom_list)


def add_ignored_atom(target_atom, atom, num=1):
    target_atom.belong.formula[atom] += 1
    for _ in range(num):
        target_atom.append_atom(atom)
    add_feature(target_atom, '-' + str(num) + atom)


def del_ignored_atom(target_atom, atom):
    target_atom.belong.formula[atom] -= 1
    del_feature(target_atom, '-' + str(target_atom.bond_list.count(atom)) + atom)
    while atom in target_atom.bond_list:
        target_atom.remove_atom(atom)
        # if single:
        #     break


'''
def inquire_id(atom_feature: int):
    return (atom_feature & (((1 << (MAX_FEATURE_NUM + MAX_ID_NUM)) - 1) & (~((1 << MAX_FEATURE_NUM) - 1))))>>MAX_FEATURE_NUM
# 用于实现对id的查询,至于为什么不用字符串切片,大概是因为位运算它快吧
# 此函数已正式弃用,但因为位运算很帅,所以将其保留
'''


def inquire_feature(atom_feature: int, feature: str):
    return atom_feature & (1 << HASH_FEATURE_TABLE[feature])


# def inquire_features(atom_feature: int, feature: str):
#     return (atom_feature & (((1 << (HASH_FEATURE_TABLE[feature] + 3)) - 1) & (~((1 << HASH_FEATURE_TABLE[feature]) - 1)))) >> HASH_FEATURE_TABLE[feature]
#

def add_feature(atom: Atom, feature: str):
    atom.feature |= (1 << HASH_FEATURE_TABLE[feature])


def del_feature(atom: Atom, feature: str):
    atom.feature &= ~(1 << HASH_FEATURE_TABLE[feature])


def del_atom(self):
    del_target = self
    visit_point = 0
    will_visit = [self]
    while visit_point < len(will_visit):
        self = will_visit[visit_point]
        if self.bond_list.count(del_target) >= 1:
            feature_name = "-" + str(self.bond_list.count(del_target)) + del_target.name
            del_feature(self, feature_name)
            while self.bond_list.count(del_target) >= 1:
                self.bond_list[self.bond_list.index(del_target)] = '-1h'
        for connect_index in range(len(self.bond_list)):
            if (self.bond_list[connect_index] == '-1h') or (self.bond_list[connect_index] in will_visit):
                continue
            elif type(self.bond_list[connect_index]) is list:
                if del_target in self.bond_list[connect_index]:
                    if del_target.name == 'c':
                        del_feature(self, 'c6-')
                    else:
                        del_feature(self, 'no2-')
                    self.bond_list[connect_index] = '-1h'
            else:
                will_visit.append(self.bond_list[connect_index])

        visit_point += 1
    del del_target


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
            if (isinstance(fa, str)) or (fa in will_visit) or (type(fa) is list):
                continue
            else:
                will_visit.append(fa)
                ansL.append(str(distant) + str(fa.feature))
        visit_point += 1
    ansL.sort()
    return '.'.join(ansL)


def connect(target_atom_list, is_cyclization=False):
    for tmp in range(len(target_atom_list) - 1):
        add_bond(target_atom_list[tmp], target_atom_list[tmp + 1])
    if is_cyclization:
        add_bond(target_atom_list[0], target_atom_list[-1])


if __name__ == '__main__':
    print(MAX_FEATURE_NUM, HASH_FEATURE_TABLE)
    # a = Molecule()
    # b = Molecule()
    # c1 = Atom('c', a)
    # c2 = Atom('c', a)
    # c3 = Atom('c', b)
    # c4 = Atom('c', b)
    # c1.add_bond(c2, 2)
    # c3.add_bond(c4)
    # c3.add_bond(c4)
    # a.update()
    # b.update()
    # print(a.feature)
    # print(b.feature)
    # print(a == b)

    benzaldehyde = Molecule()
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

    connect([c2, c3, c4, c5, c6, c7], True)
    add_bond(c1, c4)
    add_bond(c1, o1, 2)
    # c1.add_bond(o1)
    # c1.add_bond(o1)

    benzaldehyde1 = Molecule()

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

    connect([c11, c21, c31, c41, c51, c611, c71])
    add_bond(c11, c611)
    add_bond(c71, o11, 2)

    benzaldehyde.update()
    benzaldehyde1.update()
    print(benzaldehyde.feature)
    print(benzaldehyde1.feature)
    print()
    print(bool(benzaldehyde.feature == benzaldehyde1.feature))
    for i in c1.bond_list:
        if i == '-1h':
            continue
        print(i.overall_f)

    add_pi_bond([c2, c3, c4, c5, c6, c7])
    print(c1.name, c1.overall_f)
    print()
    print(benzaldehyde.composition)
    for i in benzaldehyde.composition:
        print(i.name, i.feature, i.overall_f)

    print(c1.bond_list, c1)
    print(c2.bond_list, c2.name)
    print()
    print(benzaldehyde.unsaturation)
    # data = json.dumps(benzaldehyde.feature)
    # print(data)
