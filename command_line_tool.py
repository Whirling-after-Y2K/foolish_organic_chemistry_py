import json

import organic_chemistry as oc

SAVE_PATH = "molecules"


class CommandLine:
    def __init__(self):
        self.molecule_dict = {}
        self.current_molecule = None
        self.current_name = None

    def new_molecule(self, name):
        if name in self.molecule_dict:
            print("此名称已被占用")
        else:
            self.molecule_dict[name] = [oc.Molecule()]
            self.current_name = name
            self.current_molecule = self.molecule_dict[name][0]

    def change_molecule(self, name=''):
        if name == '':
            print(self.current_name)
        else:
            self.current_molecule = self.molecule_dict[name][0]
            self.current_name = name

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
            print(f'index:{atom_index}')
            print(f'{atom.name}:')
            for connect in atom.bond_list:
                if connect is None:
                    break
                print(f"\t-{connect.name}")

    def save(self):
        self.current_molecule.update()

    def save_to_local(self):
        data = json.dumps(self.current_molecule.feature)
        with open(SAVE_PATH + '\\' + self.current_name + '.json', 'w') as f:
            f.write(data)


# print(oc.data)
a = CommandLine()
a.new_molecule("1")
a.add_atom('c',4)
a.save()
a.save_to_local()