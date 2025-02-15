import json
import sys

import organic_chemistry as oc

SAVE_PATH = "molecules"


class Console:
    def __init__(self):
        self.molecule_dict = {}
        self.current_molecule = None
        self.current_name = None
        self.log = ''

    def new_molecule(self, name):
        if name in self.molecule_dict:
            print("此名称已被占用")
        else:
            self.molecule_dict[name] = [oc.Molecule()]
            self.current_name = name
            self.current_molecule = self.molecule_dict[name][0]

    def change_molecule(self, name=''):
        if name in self.molecule_dict:
            self.current_molecule = self.molecule_dict[name][0]
            self.current_name = name
        else:
            print("此名称不存在")

    def rename(self, new_name):
        if new_name in self.molecule_dict:
            print("此名称已被占用")
        else:
            original_name = self.current_name
            self.molecule_dict[new_name] = self.molecule_dict[original_name]
            self.current_name = new_name
            del self.molecule_dict[original_name]

    def is_equal(self, name1, name2):
        return self.molecule_dict[name1][0].feature == self.molecule_dict[name2][0].feature

    def add_atom(self, atom_name, num=1):
        for _ in range(num):
            self.molecule_dict[self.current_name].append(oc.Atom(atom_name, self.current_molecule))

    def del_atom(self, atom_index):
        oc.del_atom(self.molecule_dict[self.current_name][atom_index])
        # print(self.molecule_dict[self.current_name])
        del self.molecule_dict[self.current_name][atom_index]
        # print(self.molecule_dict[self.current_name])

    def list_molecule(self):
        for i in self.molecule_dict:
            print(i)

    def list_atoms(self):
        for atom_index in range(1, len(self.molecule_dict[self.current_name])):
            atom = self.molecule_dict[self.current_name][atom_index]
            print(f'index:{atom_index}', end='\t')
            print(f'{atom.name}:')
            for connect in atom.bond_list:
                if type(connect) is str:
                    print('\t' + connect)
                elif type(connect) is list:
                    # print(connect)
                    print('\tpi:',end='\n\t')
                    for i in connect:
                        print(f"\t-{i.name} index:{self.molecule_dict[self.current_name].index(i)}", end=' ')
                    print()
                else:
                    print(f"\t-{connect.name} index:{self.molecule_dict[self.current_name].index(connect)}")

    def connect_atom(self, atom_index_list, is_cyclization=False, num=1):
        if len(atom_index_list) == 2:
            atom1 = self.molecule_dict[self.current_name][atom_index_list[0]]
            atom2 = self.molecule_dict[self.current_name][atom_index_list[1]]
            oc.add_bond(atom1, atom2, num)
        else:
            atom_list = [self.molecule_dict[self.current_name][atom_index] for atom_index in atom_index_list]
            oc.connect(atom_list, is_cyclization)

    def break_atom(self, atom_index1, atom_index2):
        # if len(atom_index_list) == 2:
        atom1 = self.molecule_dict[self.current_name][atom_index1]
        atom2 = self.molecule_dict[self.current_name][atom_index2]
        oc.break_bond(atom1, atom2)

    def connect_pi(self, atom_index_list):
        # atom1 = self.molecule_dict[self.current_name][atom_index_list[0]]
        atom_list = [self.molecule_dict[self.current_name][atom_index] for atom_index in atom_index_list]
        oc.add_pi_bond(atom_list)

    def add_ignored(self, atom_index, ignored_atom):
        oc.add_ignored_atom(self.molecule_dict[self.current_name][atom_index], ignored_atom)

    def del_ignored(self, atom_index, ignored_atom):
        oc.del_ignored_atom(self.molecule_dict[self.current_name][atom_index], ignored_atom)

    def load(self, file):
        if file.split('.')[-1] != 'txt':
            print('Unsupported file formats')

        local = open(file, 'r', encoding='utf-8')
        sys.stdin = local
        if self.current_name is not None:
            input()
        input()
        # print(sys.stdin)

    def save(self):
        self.current_molecule.update()

    def save_to_local(self):
        self.save()
        data = json.dumps(self.current_molecule.feature)
        with open(SAVE_PATH + '\\' + self.current_name + '.txt', 'w', encoding='utf-8') as f:
            f.write(data)
            f.write('\n')
            f.write(f'nm {self.current_name}\n')
            for atom_index in range(1, len(self.molecule_dict[self.current_name])):
                atom = self.molecule_dict[self.current_name][atom_index]
                f.write(f'a {atom.name}\n')
            is_visited = []
            for atom_index in range(1, len(self.molecule_dict[self.current_name])):
                atom = self.molecule_dict[self.current_name][atom_index]
                is_visited.append(atom_index)
                pi = False
                for connect in atom.bond_list:
                    if type(connect) is str:
                        if connect != '-1h':
                            f.write(f'c {atom_index} {connect}\n')
                    elif type(connect) is list:
                        # is_visited.append(connect)
                        pi = connect
                        continue
                    elif self.molecule_dict[self.current_name].index(connect) in is_visited:
                        continue
                    else:
                        f.write(f'c {atom_index},{self.molecule_dict[self.current_name].index(connect)}\n')
                if pi and pi not in is_visited:
                    tmp_pi = [str(self.molecule_dict[self.current_name].index(i)) for i in pi]
                    # print(','.join(tmp_pi))
                    f.write('c_pi ' + ','.join(tmp_pi) + '\n')
                    is_visited.append(pi)
            f.write('save\n')
            f.write('end\n')

    def run(self, file=''):
        console_in = sys.stdin
        if file != '':
            f = open(file, 'r', encoding='gbk')
            sys.stdin = f
        while True:
            user_input = input(f'{self.current_name}>').split(' ')
            command = user_input[0]

            tmp = 0
            while tmp < len(user_input):
                if user_input[tmp]:
                    tmp += 1
                else:
                    user_input.pop(tmp)
            self.log += ' '.join(user_input)
            self.log += '\n'
            # print(f'your input:{user_input}')

            try:
                match command:
                    case 'q' | 'quit':
                        break
                    case 's' | 'save':
                        self.save()
                    case 'save_to_local' | 'stl':
                        self.save_to_local()
                    case 'c' | 'ca' | 'connect_atom':
                        atom_index_list = [int(atom_index) for atom_index in user_input[1].split(',')]
                        if len(user_input) == 3:
                            if user_input[2].isdigit():
                                self.connect_atom(atom_index_list, False, int(user_input[2]))
                            elif user_input[2] in ['t', 'true', 'True']:
                                self.connect_atom(atom_index_list, True)
                            else:
                                self.add_ignored(atom_index_list[0], user_input[2])
                        else:
                            self.connect_atom(atom_index_list)
                    case 'b':
                        atom_index_list = [int(atom_index) for atom_index in user_input[1].split(',')]
                        self.break_atom(int(atom_index_list[0]), int(atom_index_list[1]))
                    case 'c_pi' | 'connect_pi':
                        atom_index_list = [int(atom_index) for atom_index in user_input[1].split(',')]
                        self.connect_pi(atom_index_list)
                    case 'cm' | 'change_molecule':
                        self.change_molecule(user_input[1])
                    case 'print_name' | 'current':
                        print(self.current_name)
                    case 'la' | 'list' | 'list_atoms':
                        self.list_atoms()
                    case 'lm' | 'list_molecule':
                        self.list_molecule()
                    case 'rn' | 'rename':
                        self.rename(user_input[1])
                    case 'nm' | 'new_molecule':
                        self.new_molecule(user_input[1])

                    case 'a' | 'add' | 'add_atom':
                        if len(user_input) == 3:
                            self.add_atom(user_input[1], int(user_input[2]))
                        else:
                            self.add_atom(user_input[1])
                    case 'd' | 'del':
                        if user_input[-1].isdigit():
                            self.del_atom(int(user_input[-1]))
                        else:
                            self.del_ignored(int(user_input[1]), user_input[2])
                    case 'save_log':
                        with open('log.txt', 'w') as f:
                            f.write(self.log)
                            f.write('q\n')
                    case 'end_of_load':
                        sys.stdin = original_in
                        # f.close()
                    case 'end':
                        sys.stdin = console_in
                    case 'load':
                        original_in = sys.stdin
                        self.load(user_input[-1])
                        # sys.stdin = orginal_in
                    case '==' | 'is_equal':
                        if self.is_equal(user_input[1], user_input[2]):
                            print('True')
                        else:
                            print('False')
                    case _:
                        print('wrong command')
            except IndexError:
                print('wrong command')
                raise


if __name__ == '__main__':
    # print(oc.data)
    a = Console()
    a.run('log.txt')
    # a.run()
