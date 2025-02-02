import json

CHEMISTRY_BOND_DICT = {'c': 4, 'o': 2, 'n': 3}
HASH_FEATURE_TABLE = {}
SAVE_NAME_NUM = 4
MAX_FEATURE_NUM = SAVE_NAME_NUM  # 后4位用于储存self_name 元素名 最多16个

with open("chemistry_feature.json", 'r', encoding='utf-8') as f:
    raw_feature_table = json.loads(f.read())
del f

for i in raw_feature_table:
    if i == "comment":
        continue
    elif i == "self_name":
        tmp_index2 = 0
        for j in raw_feature_table[i]:
            HASH_FEATURE_TABLE[j[0]] = tmp_index2  # tmp_index2用于为每一种名称编号 表示储存在后4位中该元素为tmp_index2
            HASH_FEATURE_TABLE[tmp_index2] = j[0]
            tmp_index2 += 1
        del tmp_index2
    else:
        tmp_index = SAVE_NAME_NUM
        tmp_feature_num = SAVE_NAME_NUM
        for j in raw_feature_table[i]:
            if '1' in j[0]:
                HASH_FEATURE_TABLE[j[0]] = tmp_index  # 对于只占一个键位的特征 用三位进行储存 储存在倒数第tmp_index+1至倒数第tmp_index+3位
                tmp_feature_num += 3
                tmp_index += 3
            else:
                HASH_FEATURE_TABLE[j[0]] = tmp_index  # tmp_index用于为每一种特征编号 表示其储存在倒数第tmp_index+1位
                tmp_feature_num += 1
                tmp_index += 1

        # 每多一种特征就需要多特定位进行储存,但不同元素的特征不可能共存 tmp_feature_num用于比较至少需要多少位进行储存
        MAX_FEATURE_NUM = max(tmp_feature_num - 1, MAX_FEATURE_NUM)

del tmp_index
del raw_feature_table


class Molecule:
    def __init__(self):
        # self.name = name
        self.composition = []
        self.feature = []
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
        self.bond_list = [None for _ in range(CHEMISTRY_BOND_DICT[name])]
        self.feature = 1 << MAX_FEATURE_NUM
        self.feature |= HASH_FEATURE_TABLE[self.name]
        self.overall_f = ''
        self.belong = molecule
        molecule.composition.append(self)

    def add_bond(self, target_atom, connect_num=1, is_first=True):
        if target_atom in self.bond_list:
            feature_name = self.name + "-" + str(self.bond_list.count(target_atom)) + target_atom.name
            del_feature(self, feature_name)
            feature_name = self.name + "-" + str(self.bond_list.count(target_atom) + connect_num) + target_atom.name
        else:
            feature_name = self.name + "-" + str(connect_num) + target_atom.name
        add_feature(self, feature_name)
        if is_first:
            target_atom.add_bond(self, connect_num, False)

        add_point = self.bond_list.index(None)
        # print(add_point)
        for _ in range(connect_num):
            self.bond_list[add_point] = target_atom
            add_point += 1
            # print(add_point)

    def add_pi_bond(self, target_atom_list: list, is_first=True):  # 注意target_element_list包含自身
        tmp_list: list = target_atom_list.copy()
        self.bond_list[self.bond_list.index(None)] = tmp_list
        f = True
        for i in target_atom_list:
            if i.name == 'c':
                continue
            else:
                f = False
                add_feature(self, 'no2-')
        if f:
            add_feature(self, 'c6-')
        if is_first:
            target_atom_list.remove(self)
            for target_atom in target_atom_list:
                target_atom.add_pi_bond(tmp_list, False)


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
    if '1' in feature:
        atom.feature += (1 << HASH_FEATURE_TABLE[feature])
    else:
        atom.feature |= (1 << HASH_FEATURE_TABLE[feature])


def del_feature(atom: Atom, feature: str):
    if '1' in feature:
        atom.feature -= (1 << HASH_FEATURE_TABLE[feature])
    else:
        atom.feature &= ~(1 << HASH_FEATURE_TABLE[feature])


def del_atom(self):
    del_target = self
    visit_point = 0
    will_visit = [self]
    while visit_point < len(will_visit):
        self = will_visit[visit_point]
        for connect_index in range(len(self.bond_list)):
            if (self.bond_list[connect_index] is None) or (self.bond_list[connect_index] in will_visit):
                pass
            elif type(self.bond_list[connect_index]) is list:
                if del_target in self.bond_list[connect_index]:
                    del self.bond_list[connect_index]
            else:
                if self.bond_list[connect_index] is del_target:
                    del self.bond_list[connect_index]
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
    c1.add_bond(c4)
    c1.add_bond(o1, 2)
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
    c11.add_bond(c611)
    c71.add_bond(o11, 2)

    benzaldehyde.update()
    benzaldehyde1.update()
    print(benzaldehyde.feature)
    print(benzaldehyde1.feature)
    print()
    print(bool(benzaldehyde.feature == benzaldehyde1.feature))
    for i in c1.bond_list:
        if i is None:
            continue
        print(i.overall_f)

    c2.add_pi_bond([c2, c3, c4, c5, c6, c7])
    print(c1.name, c1.overall_f)
    print()
    print(benzaldehyde.composition)
    for i in benzaldehyde.composition:
        print(i.name, i.feature, i.overall_f)

    print(c1.bond_list, c1)
    print(c2.bond_list, c2.name)
    print()
    data = json.dumps(benzaldehyde.feature)
    print(data)
