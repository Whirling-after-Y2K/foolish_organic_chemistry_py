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

    def rename(self, original_name, new_name):
        if new_name in self.molecule_dict:
            print("此名称已被占用")
        else:
            self.molecule_dict[new_name] = self.molecule_dict[original_name]
            del self.molecule_dict[original_name]

    def add_atom(self, atom_name, num=1):
        for _ in range(num):
            self.molecule_dict[self.current_name].append(oc.Atom(atom_name, self.current_molecule))

    def list_atoms(self):
        for atom_index in range(1, len(self.molecule_dict[self.current_name])):
            atom = self.molecule_dict[self.current_name][atom_index]
            print(f'index:{atom_index}', end='\t')
            print(f'{atom.name}:')
            for connect in atom.bond_list:
                if connect is None:
                    break
                print(f"\t-{connect.name} index:{self.molecule_dict[self.current_name].index(connect)}")

    def connect_atom(self, atom_index_list, is_cyclization=False, num=1):
        if len(atom_index_list) == 2:
            atom1 = self.molecule_dict[self.current_name][atom_index_list[0]]
            atom2 = self.molecule_dict[self.current_name][atom_index_list[1]]
            atom1.add_bond(atom2, num)
        else:
            atom_list = [self.molecule_dict[self.current_name][atom_index] for atom_index in atom_index_list]
            oc.connect(atom_list, is_cyclization)

    def connect_pi(self, atom_index_list):
        atom1 = self.molecule_dict[self.current_name][atom_index_list[0]]
        atom_list = [self.molecule_dict[self.current_name][atom_index] for atom_index in atom_index_list]
        atom1.add_pi_bond(atom_list)

    def save(self):
        self.current_molecule.update()

    def save_to_local(self):
        self.save()
        data = json.dumps(self.current_molecule.feature)
        with open(SAVE_PATH + '\\' + self.current_name + '.json', 'w') as f:
            f.write(data)

    def load(self, file):
        if file.split('.')[-1] != 'json':
            print('Unsupported file formats')
        with open(file, 'r') as f:
            data = json.load(f)
            data = data[0]
            for i in data.split('.'):
                feature = bin(int(i[1::]))[2::]
                element_name = oc.HASH_FEATURE_TABLE[feature[-4::]]

    def run(self, file=''):
        if file != '':
            f = open(file, 'r')
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
                                self.connect_atom(atom_index_list)
                        else:
                            self.connect_atom(atom_index_list)
                    case 'c_pi' | 'connect_pi':
                        atom_index_list = [int(atom_index) for atom_index in user_input[1].split(',')]
                        self.connect_pi(atom_index_list)
                    case 'cm' | 'change_molecule':
                        self.change_molecule(user_input[1])
                    case 'print_name' | 'current':
                        print(self.current_name)
                    case 'la' | 'list' | 'list_atoms':
                        self.list_atoms()
                    case 'rn' | 'rename':
                        if len(user_input) == 3:
                            self.rename(user_input[1], user_input[2])
                        else:
                            self.rename(self.current_name, user_input[1])
                    case 'nm' | 'new_molecule':
                        self.new_molecule(user_input[1])

                    case 'a' | 'add' | 'add_atom':
                        if len(user_input) == 3:
                            self.add_atom(user_input[1], int(user_input[2]))
                        else:
                            self.add_atom(user_input[1])

                    case 'save_log':
                        with open('log.txt', 'w') as f:
                            f.write(self.log)
                            f.write('q\n')

                    case _:
                        print('wrong command')
            except IndexError:
                print('wrong command')
                raise
        if file != '':
            f.close()


if __name__ == '__main__':
    # print(oc.data)
    a = Console()
    a.run()
