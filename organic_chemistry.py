import json
from heapq import heappush

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
print(MAX_FEATURE_NUM, HASH_FEATURE_TABLE)


class Fheap:
    def __init__(self, n):
        self.heap = [None for _ in range(n)]

    def add(self, value):
        for i in range(len(self.heap)):
            if type(self.heap[i]) is Atom:
                if self.heap[i].overall_f > value.overall_f:
                    value, self.heap[i] = self.heap[i], value
            elif type(self.heap[i]) is list:
                value, self.heap[i] = self.heap[i], value
            else:
                self.heap[i] = value
                return

    def add_pi(self, value):
        self.heap[self.heap.index(None)] = value

    def inquire_none(self):
        try:
            return self.heap.index(None)
        except ValueError:
            raise

    def __iter__(self):
        return iter(self.heap)


class Molecule:
    def __init__(self, name):
        self.name = name
        self.composition = {}
        self.feature_dict = {}
        self.link_dict = {}
        self.c_num = 0
        self.o_num = 0
        self.n_num = 0

    def update_name(self, new_name):
        self.name = new_name

    def add_atom(self, atom_name, num=1):
        atom_name = atom_name.lower()

        match atom_name:
            case "c":
                atom_num = self.c_num
                self.c_num += num
            case "n":
                atom_num = self.n_num
                self.n_num += num
            case "o":
                atom_num = self.o_num
                self.o_num += num
            case _:
                raise ValueError

        for _ in range(num):
            atom = Atom(atom_name, self)
            self.composition[atom_name + str(atom_num)] = atom
            self.feature_dict[atom_name + str(atom_num)] = atom.feature
            atom_num += 1


class Atom:
    def __init__(self, name: str, molecule: Molecule):
        # global element_id
        # element_id += 1
        # if element_id > ((1 << MAX_ID_NUM) - 2):
        #     raise OverflowError
        name = name.lower()
        self.name = name
        self.bond_list = Fheap(CHEMISTRY_BOND_DICT[name])
        self.feature = 1 << MAX_FEATURE_NUM
        # self.feature += (element_id << MAX_FEATURE_NUM)
        self.feature |= HASH_FEATURE_TABLE[self.name]
        self.overall_f = str(self.feature)
        self.belong = molecule
        molecule.composition[self] = self.feature
        # self.feature = 0b00000
        # self.feature += FEATURE_TABLE[name]

    def add_bond(self, target_atom, connect_num=1, is_first=True):

        left_index = self.bond_list.inquire_none()
        feature_name = self.name + "-" + str(connect_num) + target_atom.name
        add_feature(self, feature_name)
        # update_feature(self, target_atom)

        self.belong.composition[self] = self.feature
        for _ in range(connect_num):
            self.bond_list.add(target_atom)
            left_index += 1

            # if is_first:
            #     target_atom.add_bond(self, connect_num, False)
        # self.feature += FEATURE_TABLE[str(connect_num) + target_element]

    # def __eq__(self, other):
    #     return self.feature == other.feature

    def add_pi_pond(self, target_atom_list: list, is_first=True):  # 注意target_element_list包含自身
        tmp_list: list = target_atom_list.copy()
        self.bond_list.add_pi(tmp_list)
        if is_first:
            target_atom_list.remove(self)
            for target_atom in target_atom_list:
                target_atom.add_pi_pond(tmp_list, False)


'''
def inquire_id(atom_feature: int):
    return (atom_feature & (((1 << (MAX_FEATURE_NUM + MAX_ID_NUM)) - 1) & (~((1 << MAX_FEATURE_NUM) - 1))))>>MAX_FEATURE_NUM
# 用于实现对id的查询,至于为什么不用字符串切片,大概是因为位运算它快吧
'''


# def update_id()

def inquire_feature(atom_feature: int, feature: str):
    return atom_feature & (1 << HASH_FEATURE_TABLE[feature])


def add_feature(atom: Atom, feature: str):
    atom.feature |= (1 << HASH_FEATURE_TABLE[feature])


def del_feature(atom: Atom, feature: str):
    atom.feature &= ~(1 << HASH_FEATURE_TABLE[feature])


# def update_feature(atom: Atom, target_atom: Atom):
#     is_used = []
#     for indirect_index in range(len(target_atom.bond_list)):
#         if target_atom.bond_list[indirect_index] is None:
#             break
#         else:
#             if target_atom.bond_list[indirect_index] in is_used:
#                 continue
#             target_atom.overall_f += update_feature(target_atom, target_atom.bond_list[indirect_index])
#             is_used.append(target_atom.bond_list[indirect_index])
#     return target_atom.overall_f


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
    for i in c1.bond_list:
        if i is None:
            continue
        print(i.overall_f)

    c2.add_pi_pond([c2, c3, c4, c5, c6, c7])
    print(c1.name, c1.overall_f)
    print()
    print(benzaldehyde.composition)
    for i in benzaldehyde.composition:
        print(i.name, benzaldehyde.composition[i], i)
    # for i in c1.bond_list:
    #     if type(i) is Element:
    #         print(i.name, i.feature, id(i))
    # print()
    # for i in c2.bond_list:
    #     if type(i) is Element:
    #         print(i.name, i.feature, id(i))
    #     elif type(i) is list:
    #         print()
    #         for j in i:
    #             print(j.name, j.feature, id(j))
    print(c1.bond_list, c1)
    print(c2.bond_list, c2.name)
    # print(element_id)
