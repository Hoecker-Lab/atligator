import io
import pickle
import sys
from argparse import FileType, Namespace
from tkinter import font, Tk, Button, ttk, Label, StringVar, BooleanVar, filedialog, Checkbutton, Radiobutton, IntVar, \
    Spinbox, Frame, W, messagebox, Toplevel
from typing import Dict, Any, List, Tuple

from atligator.acomplex import read_complex_descriptor, ComplexDescriptor
from atligator.argument_collection import help_txts, arg
from atligator.pdb_util import canonical_amino_acids as aa


def next_row(current: int, value: int = None):
    """
    Returns an incremented version of current value. if value is set the sum of current and value is returned.
    Is used for row numbers for grid window manager.
    :param current: Is returned. Either +1 or +value if value is stated.
    :param value: Added to current if stated
    :return: current + 1 or current + value (if stated)
    """
    if value is None:
        return 1 + current
    else:
        return current + value


def normalize_couples(couples: List[List[str]]) -> List[List[str]]:
    """
    Will normalize a coupling list. This list consists of lists with strings. These inner lists represent the rows
    in coupling section of the descriptor. Thus, all residues represented by strings in one row are coupled in mutation.
    Normalization will add empty strings ("") in these row lists so the residue numbers increase in the right order:
    Columns, then rows. For example:
    1 5          1 ''  11
    2 6         ''  6  13
    3 7   or     3 ''  15        are saved as [['1', '5'],['2', '6'],['3', '7'],['4', '8']]
    4 8          4  8  ''                  or [['1', '', '11'],['', '6', '13'],['3', '', '15'],['4', '8', '']]
    :param couples: List of list with strings representing residues or wildcards ('')
    :return: Same list with inserted wildcards to follow the right order (increasing in columns then rows)
    """

    def sort_row(list_list: List[List[str]]) -> List[List[str]]:
        """
        Sorts rows of a matrix (list of list with int as strings) and returns it with strings again
        :param list_list: list of lists to sort
        :return: sorted lists in list
        """
        for r_i, r in enumerate(list_list):
            for x in range(len(r)):
                r[x] = int(r[x])
            list_list[r_i] = sorted(r)
        for r in list_list:
            for y in range(len(r)):
                r[y] = str(r[y])
        return list_list

    def next_non_zero(list_list: List[List[str]], col: int, row_i: int, x_y: str = "x") -> int or float:
        """
        Returns the next residue id that is not zero. If the next one is zero it walks in a certain direction (x_y)
        and returns the next non zero value or an incredibly high value if it doesn't find any.
        :param list_list: List of lists to look for this residue id
        :param col: starting point: column number
        :param row_i: starting point: row number
        :param x_y: direction to walk
        :return: int value of residue id
        """
        value = 0
        j = 0
        while value == 0:
            if x_y == "x":
                value = int(list_list[row_i][col + j])
            elif x_y == "y":
                if row_i + j == len(list_list) or col >= len(list_list[row_i + j]):
                    return j, 10000
                value = int(list_list[row_i + j][col])
            else:
                sys.exit()
            j += 1
        return j - 1, value

    couples = sort_row(couples)
    if all(e == "" for e in couples):
        return [["", ""]]
    not_finished = True
    max_len = 0
    while not_finished:
        not_finished = False
        max_len = max(len(l) for l in couples)
        for u in range(max_len):
            sth_happens = True
            while sth_happens:
                sth_happens = False
                for i, row in enumerate(couples):
                    if i != len(couples) - 1 and u < len(row):
                        c, this_res = next_non_zero(couples, u, i)
                        c, next_res = next_non_zero(couples, u + c, i + 1, "y")
                        if next_res < this_res:
                            row.insert(u, "0")
                            sth_happens = True
                            if len(row) > max_len:
                                not_finished = True
                    elif u < len(row):
                        c, this_res = next_non_zero(couples, u, i)
                        c, next_res = next_non_zero(couples, u + 1 + c, 0, "y")
                        if next_res < this_res:
                            row.insert(u, "0")
                            sth_happens = True
                            if len(row) > max_len:
                                not_finished = True
    for row in couples:
        for _ in range(max_len - len(row)):
            row.append("")
        for res in range(len(row)):
            if row[res] == '0':
                row[res] = ''
    return couples


class Interface:
    def __init__(self, title, main: bool = True):
        if not main:
            self.window = Toplevel()
        else:
            self.window = Tk()
        default_font = font.nametofont("TkDefaultFont")
        default_font.configure(size=12, weight='bold')
        self.bold = font.Font(weight='bold')
        self.window.option_add("*Font", default_font)
        # self.window.configure(background='#e0ffec')
        self.window.title(title)
        # centers the self.window
        # self.window.eval('tk::Place.self.window %s center' % self.window.winfo_pathname(self.window.winfo_id()))
        # self.window.eval('tk::PlaceWindow %s center' % self.window.winfo_toplevel())
        # self.window.geometry('350x200')
        self.row_ = 0

    def run_mode(self):
        pass

    def build_main_window(self, **kwargs):
        pass

    def send_cmd(self):
        pass

    @staticmethod
    def return_default(_) -> Any:
        return None

    def exit_ui(self):
        if isinstance(self, ParameterManager):
            del self.window
            sys.exit()
        else:
            self.window.destroy()

    def insert_sendcancel(self, frame: Frame or Tk or Toplevel = None):
        if frame is None:
            frame = self.window
        stop = Button(frame, text="Cancel", command=self.exit_ui)
        stop.pack(side="right")
        send = Button(frame, text="Send", command=self.send_cmd)
        send.pack(side="right")

    def insert_widget_basics(self, name: str, frame_id: Frame):
        """
        Inserts default widgets for elements such as a help button and a label.
        :param self: ParameterManager object we're working in
        :param name: Name of the element (short version)
        :param frame_id: Frame to put widgets in
        """
        self.insert_help_button(name, frame_id)
        self.insert_element_label(name, frame_id)

    def insert_help_button(self, name: str, frame_id: Frame):
        """
        Inserts a help button for specified element (name) in frame (frame_id).
        :param self: ParameterManager object we're working in
        :param name: Name of the element (short version)
        :param frame_id: Frame to put widgets in
        """
        info_button = Button(frame_id, text="?", fg='white', bg='green', font=self.bold,
                             command=lambda: self.show_info(name))
        info_button.grid(column=0, row=self.row_, sticky=W)

    def insert_element_label(self, name: str, frame_id: Frame or Tk, ret: bool = False,
                             text: str = False, column=1, row=None, sticky=W, fg="black", pack: bool = False,
                             side: str = "right"):
        """
        Inserts a label for specified element (name) in frame (frame_id).
        Special cases (for filebox content or major frames) can be inserted by specifying
        an 'ret' which returns an element that should be assigned to that label
        (For example a dynamic file box content text)

        :param self: ParameterManager object we're working in
        :param name: Name of the element (short version)
        :param frame_id: Frame to put widgets in
        :param ret: Return an element to assign the Label to. If True, text and row needs to be specified, too.
        :param text: Text to insert in the Label. Only taken into account if 'ret' is true
        :param column: Column to (grid) pack the label in
        :param row: Row to (grid) pack the label in
        :param sticky: Adjusts the sticky side. If None, sticky is not enabled.
        :param fg: Adjusts the font color of the label. Only taken into account if 'ret' is true.
        :param pack: Switches to pack window manager instead of grid if enabled
        :param side: Side where to pack label.
        """
        row_ = self.row_ if row is None else row
        if ret is False:
            label_ = Label(frame_id, text="    " + name + ": ")
            if not pack:
                label_.grid(column=column, row=row_, sticky=W)
            else:
                label_.pack(side=side, anchor=sticky)
        else:
            element = Label(frame_id, text=text, fg=fg)
            if not pack:
                element.grid(column=column, row=row_, sticky=sticky)
            else:
                element.pack(side=side, anchor=sticky)
            return element

    def show_info(self, name: str):
        """
        Shows a messagebox with long id of 'name' and its help text
        :param name: name of parameter (e.g. atlas)
        """
        messagebox.showinfo('Default Info', "We all wanna change the world.")


class ParameterManager(Interface):
    def __init__(self, defaults):
        super().__init__("Welcome to atligator!")
        self.window.eval('tk::PlaceWindow %s center' % self.window.winfo_toplevel())

        self.conversion_dict: Dict[str, str] = {}
        for a in arg:
            if 'l_name' in arg[a]:
                self.conversion_dict[arg[a]['l_name']] = a
            else:
                self.conversion_dict[a] = a
        del self.conversion_dict['graphical_input']
        del self.conversion_dict['atpar']

        # Read default parameters from argument parser input
        self.defaults: Dict[str, Any] = {}
        for args in vars(defaults):
            if args not in ("graphical_input", "", "atpar"):
                self.defaults[self.conversion_dict[args]] = vars(defaults)[args]

        # If parameter is mandatory use current value as default:
        for defs in self.defaults:
            if arg[defs]["group"] == "None":
                arg[defs]["default"] = defs

        # Read descriptor as well as ligand length and residues from descriptor
        self.descriptor = read_complex_descriptor(self.defaults['descriptor'].name)
        self.lig_length = len(self.descriptor.ligand_info)
        self.residues: Dict[int, ttk.Combobox] = {i: None for i in range(self.lig_length)}
        self.coupling_info: List[Tuple[str, List[str]]] = []

        # Define groups and save status: if group frame was opened or not
        self.g_checked: Dict[str, Label] = {'mp': None, 'aa': None, 'lc': None,
                                            'fit': None, 'mut': None, 'sc': None,
                                            'cs': None, 'da': None, 'rr': None,
                                            'os': None, 'out': None}

        # Dictionary with major frames for groups
        self.maj_frames: Dict[str, List[Frame]] = {}
        # Dict for checkbutton for opening group frame. Contains label, button and current state
        self.chb: Dict[str, List[Label or Button or BooleanVar]] = {}
        for a in self.g_checked:
            self.maj_frames[a] = [None, None, None]
            self.chb[a] = [None, None, None]

        self.button_dict = {}
        self.spins: Dict[str, Spinbox or None] = {}
        self.checks: Dict[str, Checkbutton or None] = {}
        self.checks_val: Dict[str, BooleanVar or None] = {}
        self.radios_val: Dict[str, StringVar or str] = {}
        self.radios_choices: Dict[str, List[str]] = {}
        self.radios: Dict[str, List[Radiobutton or None]] = {}
        self.radio_frames: Dict[str, Frame] = {}
        self.entrys: Dict[str, Any or None] = {}
        self.files: Dict[str, Any or None] = {}
        self.files_curr: Dict[str, str or None] = {}
        for a in arg:
            if arg[a]['type'] in ('int', 'float'):
                self.spins[a] = None
            if a != "G" and arg[a]['type'] == 'bool':
                self.checks[a] = None
                self.checks_val[a] = None
            if arg[a]['type'] == 'str' and arg[a]['choices'] != "":
                self.radios_val[a] = ""
                self.radios_choices[a] = arg[a]['choices']
                self.radios[a] = [None] * len(arg[a]['choices'])
                self.radio_frames[a] = None
            if arg[a]['type'] == 'str' and arg[a]['choices'] == "" and a != "A":
                self.entrys[a] = None
            if arg[a]['type'] in ("FileType('r')", "FileType('w')") and a != "descriptor":
                self.files[a] = None
                self.files_curr[a] = None
        self.result = {}
        self.namespace = Namespace()
        self.bottom_label = Label(self.window, text="")
        self.bottom_label.pack(side="bottom")
        self.row_ = 0
        self.atpar = None
        for args in vars(defaults):
            if args == "atpar":
                if vars(defaults)[args] is not None:
                    self.atpar = vars(defaults)[args]
                break

    def run_mode(self):
        self.build_main_window()
        # print(self.namespace)
        return self.namespace

    def combine_results(self):

        def modify_descriptor(name: str = "descriptor"):
            def delete_old_coupling(l: List) -> List:
                l.reverse()
                splitter = 0
                for x in l:
                    if x.startswith("-"):
                        splitter += 1
                if splitter < 2:
                    l.insert(0, "-----------\n")
                else:
                    counter = 0
                    for x in l:
                        if x.startswith("-"):
                            break
                        counter += 1
                    if counter > 0:
                        del l[:counter]
                l.reverse()
                return l

            modify_lig_res: bool = not all(self.residues[res].get() == self.descriptor.ligand_info[res][2] for res in
                                           range(self.lig_length))
            modify_coupling: bool = len(self.coupling_info) != 0
            if modify_lig_res or modify_coupling:
                line_list = []
                with open(self.defaults[name].name, "r") as df:
                    for lines in df:
                        line_list.append(lines)
                if modify_lig_res:
                    for j in range(self.lig_length):
                        for char in range(len(line_list[j])):
                            if line_list[j][char:char + 3] in aa:
                                line_list[j] = line_list[j].replace(line_list[j][char:char + 3], self.residues[j].get())
                                break
                if modify_coupling:
                    line_list = delete_old_coupling(line_list)
                    for row in self.coupling_info:
                        for r, res in enumerate(row[1]):
                            row[1][r] = (3 - len(str(res))) * " " + str(res)
                        line_list.append(row[0] + " " + " ".join(row[1]) + "\n")
                with open(self.defaults[name].name, "w") as df:
                    df.write("".join(line_list))
                print("descriptor should have been changed")
                desc_file_func = FileType('r') if arg[name]['type'] == "FileType('r')" else FileType('w')
                self.result[name] = desc_file_func(self.defaults[name].name)
            else:
                self.result[name] = self.defaults[name]
                print("descriptor not changed")

        for i in self.defaults:
            if i == "descriptor":
                # If nothing changed in ligand residues or coupling - take the same as default
                modify_descriptor()

        for i in self.checks_val:  # Include stati of checkboxes to results
            self.result[i] = self.checks_val[i].get()
        for i in self.spins:  # Include all spinbox values to results and convert in int or float
            if "." in self.spins[i].get():
                self.result[i] = float(self.spins[i].get())
            else:
                self.result[i] = int(self.spins[i].get())
        for i in self.radios_val:  # Include the current stati of all radio values to results
            self.result[i] = self.radios_val[i].get()
        for i in self.entrys:  # Include all strings (from Comboboxes or Entrys) into results
            if self.entrys[i].get() == "value unset":
                self.result[i] = None
            else:
                self.result[i] = self.entrys[i].get()
        for i in self.files:  # Include all FileType(s) (from Comboboxes or Entrys or Paths) into results
            if self.files[i] is None:
                self.result[i] = None
            elif i is not "descriptor":
                file_func = FileType('r') if arg[i]['type'] == "FileType('r')" else FileType('w')
                self.result[i] = file_func(self.files[i])

    def send_cmd(self):
        if self.perfect():
            print("sending..")
            self.combine_results()
            for i in self.result:  # Create an argparser Namespace object with all arguments (long identifiers) + values
                long_id: str = list(self.conversion_dict.keys())[list(self.conversion_dict.values()).index(i)]
                setattr(self.namespace, long_id, self.result[i])
            self.window.withdraw()
            self.window.destroy()
        else:
            print("An error was detected. Check ligand residues or coupling.")

    def perfect(self) -> bool:
        if all(self.residues[x].get() in aa for x in self.residues):
            return True
        return False

    def save_params(self):
        """
        Saves all parameters to a pickle file in '*.atpar' format.
        All special formats (incompatible with pickle) are converted first and
        reconverted after importing (import_params).
        """
        def convert(lit: Dict):
            for i in lit:
                if isinstance(lit[i], BooleanVar):
                    lit[i] = bool(lit[i].get())
                elif isinstance(lit[i], StringVar):
                    lit[i] = lit[i].get()
                elif isinstance(lit[i], IntVar):
                    lit[i] = int(lit[i])
                elif isinstance(lit[i], Spinbox):
                    lit[i] = lit[i].get()
                elif isinstance(lit[i], ttk.Combobox):
                    lit[i] = lit[i].get()

        a = dict(self.defaults)
        b = dict(self.checks_val)
        c = dict(self.spins)
        d = dict(self.radios_val)
        e = dict(self.entrys)
        f = dict(self.files)
        g = dict(self.residues)
        h = self.coupling_info
        for part in a:
            if isinstance(a[part], io.TextIOWrapper):
                a[part] = "IOWrapper " + a[part].name
        data: List[Any] = [a, b, c, d, e, f, g]
        for lit_ in data:
            convert(lit_)
        data.append(h)
        par_file = filedialog.asksaveasfilename(filetypes=(("atpar files", "*.atpar"), ("all files", "*.*")))
        if not par_file:
            print("Saving process stopped. Enter valid name")
            self.bottom_label.configure(text="Saving process stopped. Enter valid name!",
                                        fg="red", bg="grey")
            return 0
        with open(par_file, 'wb') as wb:
            pickle.dump(data, wb)
        print("saved in", par_file)
        self.bottom_label.configure(text="Successfully saved to " + par_file + "!", fg="green", bg="grey")

    def import_params(self, init: str or bool = False):
        """
        Imports parameters from a '*.atpar' file (created with save_params).
        Reconverts special types from pickle (e.g. FileType)
        """
        if init:
            par_file = init
        else:
            par_file = filedialog.askopenfilename(filetypes=(("atpar files", "*.atpar"), ("all files", "*.*")))
        if not par_file:
            print("Import process stopped. Enter valid parameter filename")
            self.bottom_label.configure(text="Import process stopped. Enter valid parameter filename!",
                                        fg="red", bg="grey")
            return 0
        with open(par_file, 'rb') as rb:
            try:
                data2 = pickle.load(rb)
            except EOFError:
                print("Import process stopped. Enter valid parameter filename")
                self.bottom_label.configure(text="Import process stopped. Enter valid parameter filename!",
                                            fg="red", bg="grey")
                return 0
        self.defaults = data2[0]
        for check in self.checks_val:
            self.checks_val[check].set(data2[1][check])
        for spin in self.spins:
            self.spins[spin].delete(0, "end")
            self.spins[spin].insert(0, data2[2][spin])
        for radio in self.radios_val:
            self.radios_val[radio].set(data2[3][radio])
        for entry in self.entrys:
            self.entrys[entry].set(data2[4][entry])
        for file in self.files:
            self.files[file] = data2[5][file]
            if self.files_curr[file] is not None:
                self.files_curr[file].configure(text=data2[5][file])
        for part in self.defaults:
            val = self.defaults[part]
            if type(val) == str and val.split()[0] == "IOWrapper":
                file_func = FileType('r') if arg[part]['type'] == "FileType('r')" else FileType('w')
                self.defaults[part] = file_func(val.split()[1])
        for res in self.residues:
            self.residues[res].current(aa.index(data2[6][res]))
        self.coupling_info = data2[7]
        print("imported from", par_file)
        self.bottom_label.configure(text=f"Successfully imported from {par_file}!", fg="green")

    def update_dict(self):
        """
        Updates the checkbuttons for hiding or showing group specific frames.
        Used if a state is changed ('show' checkbutton is (de)activated).
        """
        for group in self.g_checked:
            self.button_dict[group] = self.chb[group][2]

    def show_hide(self, fr: Frame, name: str):
        """
        Toggles a frame to appear or disappear.
        :param fr: Frame to toggle
        :param name: name of the group. Used for a 'checked' label to appear or switch color
        """
        self.update_dict()
        on_off: BooleanVar = self.button_dict[name]
        self.g_checked[name].configure(text=" " * 42 + "[checked]")
        if on_off.get():
            self.g_checked[name].configure(fg='green')
            fr.pack(side="bottom", anchor="w")
            self.chb[name][1].configure(text=" Hide ")
        else:
            self.g_checked[name].configure(fg='grey')
            fr.pack_forget()
            self.chb[name][1].configure(text="Show")

    def browse_file(self, name: str, write: bool = False):
        """
        Opens a file dialog to browse for an existing file (read or write) or a new file (write).
        If write parameter is False, only existing files are allowed to choose.
        :param name: name of the corresponding parameter
        :param write: Enables non-existing files for writing
        """
        if write:
            file = filedialog.asksaveasfilename()
        else:
            file = filedialog.askopenfilename()
        if file:
            self.files_curr[name].configure(text=file, fg="black")
            self.files[name] = file

    def unset(self, name: str):
        """
        Sets the value of an FileType parameter to None. Mostly used for default setups of disabling functions.
        :param name: id of parameter (e.g. 'Ol').
        """
        self.files_curr[name].configure(text="Unset", fg="red")
        self.files[name] = None

    def adjust_coupling(self):
        """
        Opens coupling manager and refreshes coupling.
        """
        if len(self.coupling_info) == 0:
            ui = CouplingManager(descriptor=self.descriptor)
        else:
            ui = CouplingManager(descriptor=self.descriptor, coupled_res=self.coupling_info)
        self.coupling_info = ui.run_mode()

    def show_info(self, name: str):
        """
        Shows a messagebox with long id of 'name' and its help text
        :param name: name of parameter (e.g. atlas)
        """
        long_id: str = list(self.conversion_dict.keys())[list(self.conversion_dict.values()).index(name)]
        messagebox.showinfo('Info ' + long_id, help_txts[name])

    def insert_spinbox(self, name: str, frame_id: Frame, dec_dig: int, lim: Tuple = None):
        """
        Inserts a spinbox.
        :param self: ParameterManager object we're working in
        :param name: Name of the element (short version)
        :param frame_id: Frame to put widgets in
        :param dec_dig: Decimal digits that should be shown. Also alters the increment.
        :param lim: Tuple which includes minimum and maximum number for spinbox
        """
        self.insert_widget_basics(name, frame_id)
        default_num = IntVar()
        default_num.set(self.defaults[name])
        if lim is not None:
            start_ = lim[0]
            stop_ = lim[1]
        else:
            start_ = 0
            stop_ = 100
        self.spins[name] = Spinbox(frame_id, format="%." + str(dec_dig) + "f",
                                   increment=1.0 if dec_dig < 1 else 0.1,
                                   from_=start_, to=stop_, width=25,
                                   textvariable=default_num)
        self.spins[name].grid(column=2, row=self.row_)

    def insert_radio(self, name: str, frame_id: Frame):
        """
        Inserts several radio buttons. The number is based on the choices in self.radios_choices for 'name'.
        :param self: ParameterManager object we're working in
        :param name: Name of the element (short version)
        :param frame_id: Frame to put widgets in
        """
        self.insert_widget_basics(name, frame_id)
        self.radios_val[name] = StringVar()
        self.radios_val[name].set(self.defaults[name])
        self.radio_frames[name] = Frame(frame_id)
        for i in range(len(self.radios_choices[name])):
            self.radios[name][i] = Radiobutton(self.radio_frames[name],
                                               text=self.radios_choices[name][i],
                                               value=self.radios_choices[name][i],
                                               variable=self.radios_val[name])
            self.radios[name][i].grid(column=i, row=0)
        self.radio_frames[name].grid(column=2, row=self.row_, columnspan=10, sticky=W)

    def insert_checkbutton(self, name: str, frame_id: Frame):
        """
        Inserts a checkbutton.
        :param self: ParameterManager object we're working in
        :param name: Name of the element (short version)
        :param frame_id: Frame to put widgets in
        """
        self.insert_widget_basics(name, frame_id)
        self.checks_val[name] = BooleanVar()
        self.checks_val[name].set(self.defaults[name])
        self.checks[name] = Checkbutton(frame_id, text='Choose', var=self.checks_val[name])
        self.checks[name].grid(column=2, row=self.row_)

    def insert_coupling_option(self, frame_id: Frame):
        """
        Inserts button which triggers a coupling manager in a separate window.
        :param frame_id: Frame to put button in
        """
        coupling_button = Button(frame_id, text="Adjust coupling", command=self.adjust_coupling)
        coupling_button.grid(column=3, row=self.row_, sticky=W)

    def insert_combobox(self, name: str, frame_id: Frame, *args: str):
        """
        Inserts a Combobox.
        :param self: ParameterManager object we're working in
        :param name: Name of the element (short version)
        :param frame_id: Frame to put widgets in
        :param args: Additional choices for combobox value selection
        """
        self.insert_widget_basics(name, frame_id)
        default = self.defaults[name] if self.defaults[name] is not None else "value unset"
        values_ = [default]
        for x in args:
            values_.append(str(x))
        self.entrys[name] = ttk.Combobox(frame_id,
                                         text=name, values=values_, width=25)
        self.entrys[name].set(default)  # set the selected item
        self.entrys[name].grid(column=2, row=self.row_)

    def insert_filebox(self, name: str, frame_id: Frame, write_file: bool = False):
        """
        Inserts a Filebox with a browse button.
        :param self: ParameterManager object we're working in
        :param name: Name of the element (short version)
        :param frame_id: Frame to put widgets in
        :param write_file: Non-existing files are allowed, if enabled
        """
        self.insert_widget_basics(name, frame_id)
        # Change defaults to strings (from None or TextIOWrapper)
        # if self.defaults[name] is None:
        #    self.defaults[name] = "None"
        if isinstance(self.defaults[name], io.IOBase):
            self.defaults[name] = self.defaults[name].name
        if self.defaults[name] is None:
            self.files_curr[name] = self.insert_element_label(name, frame_id, True, column=2,
                                                              text="Unset", fg="red")
        else:
            self.files_curr[name] = self.insert_element_label(name, frame_id, True, column=2,
                                                              text=self.defaults[name])
        if write_file:
            button = Button(frame_id, text="Save as...",
                            command=lambda: self.browse_file(name, write_file))
        else:
            button = Button(frame_id, text="Browse",
                            command=lambda: self.browse_file(name))
        self.files[name] = self.defaults[name]  # set the default selected item
        button.grid(column=3, row=self.row_)
        if arg[name]['type'] in ("FileType('w')", "FileType('r')") and arg[name]['default'] in ("", "genmut.log"):
            unset = Button(frame_id, text="unset",
                           command=lambda: self.unset(name))
            unset.grid(column=4, row=self.row_)

    def insert_major_frame(self, name: str, description: str):
        """
        Inserts a major frame which contents a frame which can be enabled to show parameter selection.
        This double frame insertion is needed to get the arranger structure: grid(pack(grid)).
        :param self: ParameterManager object we're working in
        :param name: Name of the group (short version)
        :param description: Description of the group
        """
        self.maj_frames[name][0] = Frame(self.window)
        self.maj_frames[name][0].pack(side="top", fill="x")
        self.maj_frames[name][2] = Frame(self.maj_frames[name][0])
        self.maj_frames[name][2].pack(side="top", anchor="w", fill="both", expand="true")
        self.chb[name][0] = self.insert_element_label(name, self.maj_frames[name][2], ret=True, sticky="w",
                                                      text=description, side="left", pack=True)
        self.chb[name][2] = BooleanVar()
        self.chb[name][2].set(False)  # set check state
        self.maj_frames[name][1] = Frame(self.maj_frames[name][0])
        self.chb[name][1] = Checkbutton(self.maj_frames[name][2], text='Show', var=self.chb[name][2],
                                        command=lambda: self.show_hide(self.maj_frames[name][1], name),
                                        indicatoron="false")
        self.chb[name][1].pack(side="right", anchor="w")
        self.g_checked[name] = self.insert_element_label(name, self.maj_frames[name][2], ret=True, text=" " * 61,
                                                         sticky="e", pack=True, side="right")

    def insert_group(self, id_: str):
        """
        Adds all elements of a parameter group. Widget type is chosen by type key value in 'arg' dict.
        :param self: ParameterManager object we're working in
        :param id_: Name of the group (short version)
        """
        self.row_ = 0
        frame = self.maj_frames[id_][1]
        for ar in arg:
            if arg[ar]['group'] == id_:
                if arg[ar]['type'] in ('float', 'int'):
                    if 'limit' in arg[ar]:
                        self.insert_spinbox(ar, frame, 0 if arg[ar]['type'] == 'int' else 2, arg[ar]['limit'])
                    else:
                        self.insert_spinbox(ar, frame, 0 if arg[ar]['type'] == 'int' else 2)
                    self.row_ = next_row(self.row_)
                elif ar != "G" and arg[ar]['type'] == 'bool':
                    self.insert_checkbutton(ar, frame)
                    if ar == "Mu":
                        self.insert_coupling_option(frame)
                    self.row_ = next_row(self.row_)
                elif arg[ar]['type'] == 'str' and arg[ar]['choices'] != "":
                    self.insert_radio(ar, frame)
                    self.row_ = next_row(self.row_)
                elif arg[ar]['type'] == 'str':
                    if "alter" in arg[ar]:
                        self.insert_combobox(ar, frame, arg[ar]['alter'])
                    else:
                        self.insert_combobox(ar, frame)
                    self.row_ = next_row(self.row_)
                elif arg[ar]['type'] == "FileType('r')":
                    self.insert_filebox(ar, frame)
                    self.row_ = next_row(self.row_)
                elif arg[ar]['type'] == "FileType('w')":
                    self.insert_filebox(ar, frame, True)
                    self.row_ = next_row(self.row_)

    def build_main_window(self):
        """
        This method sets up all the different Elements/Widgets in group specific frames.
        The frames can be en- as well as disabled by a checkbox.
        All Elements contain a Help box to display the help text and a Label for the short name.
        """

        # Mandatory parameters (atlas, scaffold)
        self.insert_major_frame('mp', "Mandatory parameters: ")
        fr_mp = self.maj_frames['mp'][1]
        self.row_ = 0
        for a in ("atlas", "scaffold"):
            self.insert_filebox(a, fr_mp)
            self.row_ = next_row(self.row_)

        # Ligand residues (modifies descriptor!)
        self.insert_major_frame('aa', "Peptide residue selection: ")
        fr_aa = self.maj_frames['aa'][1]
        for res in range(self.lig_length):
            self.residues[res] = ttk.Combobox(fr_aa, width=4, state="readonly")
            self.residues[res]['values'] = aa
            self.residues[res].current(aa.index(self.descriptor.ligand_info[res][2]))
            self.residues[res].grid(column=res, row=0)

        # Life cycle parameters
        self.insert_major_frame('lc', "Life cycle parameters: ")
        self.insert_group('lc')

        # Fitness function parameters
        self.insert_major_frame('fit', "Fitness function parameters: ")
        self.insert_group('fit')

        # Mutation parameters
        self.insert_major_frame('mut', "Mutation parameters: ")
        self.insert_group('mut')

        # Pmatrix scoring parameters
        self.insert_major_frame('sc', "Pmatrix scoring parameters: ")
        self.insert_group('sc')

        # Clash scoring parameters
        self.insert_major_frame('cs', "Clash scoring parameters: ")
        self.insert_group('cs')

        # Deep atlas parameters
        self.insert_major_frame('da', "Deep atlas parameters: ")
        self.insert_group('da')

        # Rosetta parameters
        self.insert_major_frame('rr', "Rosetta parameters: ")
        self.insert_group('rr')

        # Osprey parameters
        self.insert_major_frame('os', "Osprey parameters: ")
        self.insert_group('os')

        # Output options
        self.insert_major_frame('out', "output options: ")
        self.insert_group('out')

        save_ = Button(self.window, text="Save parameters", command=self.save_params)
        save_.pack(side="left")
        import_ = Button(self.window, text="Import parameters", command=self.import_params)
        import_.pack(side="left")

        self.insert_sendcancel()

        if self.atpar:
            self.import_params(self.atpar)
            self.atpar = False

        self.window.mainloop()


class CouplingManager(Interface):
    """
    An interface type class that defines a coupling manager window.
    This windows allows the user to define coupling individually and send it to the genetic mutator (by the parameter
    manager).
    """

    def __init__(self, descriptor: ComplexDescriptor, coupled_res=None):
        super().__init__("Mutation coupling manager", main=False)
        self.coupler: Dict[int, Dict[int, Any]] = {}
        self.descriptor = descriptor
        if coupled_res is None:  # Coupling manager is opened for first time -> coupled_res is read from descriptor
            self.coupled_res: List[List[str]] = self.read_descriptor("coupled_res")
        else:  # Coupling manager was open before -> coupled_res is read from coupled_res parameter
            cs: List[List[str]] = []
            for c in coupled_res:
                cs.append(list(filter(lambda a: a != "", c[1])))
            self.coupled_res = normalize_couples(cs)
        self.binder_res: List[str] = self.read_descriptor("binder_res")
        self.max_couplers = max(len(l) for l in self.coupled_res)
        self.coupling_info = self.descriptor.coupling_info
        self.cp_frame = None
        self.cp_daughter: Dict[int, Frame] = {}

    def run_mode(self) -> List[Tuple[str, List[str]]]:
        """
        builds the main window (coupling manager window) and returns coupling info, if main window is destroyed.
        :return: couling information which was defined in read_current_coupling
        """
        self.build_main_window()
        return self.coupling_info

    def read_descriptor(self, mode: str) -> List[List[str]] or List[str]:
        """
        Reads descriptor and saves binder or coupling info depending on mode (either "coupled_res" or "binder_res").
        :param mode: decides to return binder or coupling info
        :return: descriptor info depending on mode
        """
        if mode == "coupled_res":
            couples: List[List[str]] = []
            for couple in self.descriptor.coupling_info:
                couples.append(couple[1])

            couples = normalize_couples(couples)
            return couples
        elif mode == "binder_res":
            binder_res: List[str] = []
            for binder in self.descriptor.binder_info:
                binder_res.append(binder[1])
            binder_res.insert(0, "")
            return binder_res

    def read_current_coupling(self) -> bool:
        """
        Reads current status of Comboboxes that define coupling and combines it to self.coupling_info.
        If one residue number is found twice False is returned, otherwise True.
        :return: False if one residue id is taken twice, else True.
        """

        def remove_empty_lines(cci: List[List[str or List[str]]]) -> List[List[str or List[str]]]:
            """
            Removes empty rows or columns (filled with "") when sending coupling information.
            :param cci: current coupling info which will be used by ParameterManager
            :return: cci list without empty lines
            """
            line_to_remove = []
            for line in range(len(cci)):
                if all(element == "" for element in cci[line][1]):
                    line_to_remove.append(line)
            for ltr in sorted(line_to_remove, reverse=True):
                del cci[ltr]
            col_to_remove = []
            for col in range(len(cci[0][1])):
                if all(cci[element][1][col] == "" for element in range(len(cci))):
                    col_to_remove.append(col)
            for ctr in sorted(col_to_remove, reverse=True):
                for line in range(len(cci)):
                    del cci[line][1][ctr]
            return cci

        doublecheck = set()
        current_coupling_info: List[List[str, List[str]]] = [None for _ in range(len(self.coupler))]
        for i in self.coupler:
            s = [None for _ in range(len(self.coupler[i]))]
            if len(self.descriptor.coupling_info) > 0:
                # Only chain of first couple is used
                row: List[str or List[str]] = [self.descriptor.coupling_info[0][0], s]
            else:  # If no couples are there use chain of first binder
                row: List[str or List[str]] = [self.descriptor.binder_info[0][0], s]
            for j in self.coupler[i]:
                row[1][j] = self.coupler[i][j].get()
                if row[1][j] not in doublecheck:
                    doublecheck.add(row[1][j])
                elif row[1][j] != "":
                    print("One residues is coupled twice!")
                    messagebox.showinfo('Sending not possible', "One residues is coupled twice!")
                    return False
            current_coupling_info[i] = row
        current_coupling_info = remove_empty_lines(current_coupling_info)
        for i in range(len(current_coupling_info)):
            current_coupling_info[i] = tuple(current_coupling_info[i])
        self.coupling_info = current_coupling_info
        return True

    def send_cmd(self):
        """
        Reads current coupling status and destroys the window.
        Is used when "send" button is hit to save coupling into the ParameterManager
        """
        print("Refreshing coupling data..")
        if self.read_current_coupling():
            self.window.destroy()

    def add_del_cr(self, name: str):
        """
        Adds or deletes a row or column in coupling matrix.
        :param name: defines the behaviour. There are four possibilities: "+col", "-col", "+row" or "-row".
        :return:
        """

        def ad_row(_self, add_or_del: str):
            """
            For adding or deleting a row.
            :param _self: self CouplingManager object instance
            :param add_or_del: "+" or "-"
            """
            if add_or_del == "+":
                m = max(_self.coupler.keys())
                d = [None for _ in _self.coupled_res[m]]
                for i, j in enumerate(_self.coupled_res[m]):
                    d[i] = ""
                _self.coupled_res.append(d)
                _self.load_couples(_self.cp_frame, m + 1)

            elif add_or_del == "-":
                m = max(_self.coupler.keys())
                if m >= 1:
                    _self.cp_daughter[m].destroy()
                    del _self.coupler[m]
                    del _self.cp_daughter[m]
                    del _self.coupled_res[m]
                # del _self.coupler[m]

        def ad_col(_self, add_or_del: str):
            """
            For adding or deleting a column.
            :param _self: self CouplingManager object instance
            :param add_or_del: "+" or "-"
            """
            if add_or_del == "+":
                bb = _self.binder_res
                m = max(_self.coupler[0].keys()) + 1  # works if at least one line is there
                for i in range(len(_self.coupled_res)):
                    _self.coupler[i][m] = ttk.Combobox(_self.cp_daughter[i], width=4, state="readonly")
                    _self.coupler[i][m]['values'] = bb
                    _self.coupler[i][m].current(0)
                    _self.coupler[i][m].grid(column=m + 1, row=0)
                    _self.coupled_res[i].append("")
                _self.max_couplers += 1

            elif add_or_del == "-":
                m = max(_self.coupler[0].keys())
                if m > 1:
                    for i in range(len(_self.coupled_res)):
                        _self.coupled_res[i].pop()
                        _self.coupler[i][m].destroy()
                        del _self.coupler[i][m]
                    _self.max_couplers -= 1

        pm = name[0]
        cr = name[1:]
        if cr == "row":
            ad_row(self, pm)
        elif cr == "col":
            ad_col(self, pm)
        else:
            sys.exit()

    def insert_plusminus(self, frame_id: Frame, row_or_column: str, plus_or_minus: str):
        """
        Inserts a button in desired frame. This button either adds or deletes a row or column in coupling matrix.
        The behaviour is defined in row_or_column or plus_or_minus.
        :param frame_id: The frame it should be placed in
        :param row_or_column: Defines the axis the button should work on ("col" columns or "row" rows in matrix)
        :param plus_or_minus: Defines the action. Either "+" for adding or "-" for deleting a line.
        """
        if row_or_column == "col":
            side = "right"
            text = " " + plus_or_minus + " "
        elif row_or_column == "row":
            side = "bottom"
            text = plus_or_minus
        else:
            sys.exit()
        info_button = Button(frame_id, text=text, fg='white', bg='bisque', font=self.bold,
                             command=lambda: self.add_del_cr(plus_or_minus + row_or_column))
        info_button.pack(side=side)

    def load_couples(self, frame: Frame or Toplevel, single: int = None):
        """
        Loads coupling matrix which harbours rows and columns that define residues groups which shall be coupled.
        By default all residues in self.coupled_res (defined from descriptor in __init__()) are loaded and
        build in the window.
        If single is stated, only this single row position is reloaded (e.g. if a new row is created with '+' button)
        :param frame: Frame the matrix should be placed in
        :param single: Defines single row to (re)build
        """
        bb = self.binder_res
        if single is None:
            rows = range(len(self.coupled_res))
        else:
            rows = range(single, single + 1)
        for i in rows:
            self.cp_daughter[i] = Frame(frame)
            self.cp_daughter[i].pack(side="top")
            num = str(i + 1) if len(str(i + 1)) >= 2 else "0" + str(i + 1)
            self.insert_element_label(f"Couple #{num}", self.cp_daughter[i], row=0, column=0)
            self.coupler[i] = {}
            for j in range(self.max_couplers):
                self.coupler[i][j] = ttk.Combobox(self.cp_daughter[i], width=4, state="readonly")
                self.coupler[i][j]['values'] = bb
                self.coupler[i][j].current(bb.index(self.coupled_res[i][j]))
                self.coupler[i][j].grid(column=j + 1, row=0)

    def build_main_window(self):
        """
        Builds the main window - in this case a coupling managing window.
        The window will last and block commands after this one until it's destroyed.
        """
        main_top = Frame(self.window)
        main_top.pack(side="top", fill="x", expand=False)
        main_bottom = Frame(self.window)
        main_bottom.pack(side="bottom", fill="both", expand=True)
        main_bottom_top = Frame(main_bottom)
        main_bottom_top.pack(side="top", fill="both", expand=True)

        main_left = Frame(main_bottom_top)
        main_left.pack(side="left", fill="y", expand=False)
        main_right = Frame(main_bottom_top)
        main_right.pack(side="right", fill="both", expand=True)

        cols_left = Frame(main_top)
        cols_left.pack(side="left", expand=False)
        cols_right = Frame(main_top)
        cols_right.pack(side="right", fill="x", expand=True)
        self.insert_help_button("help", cols_left)
        for i, pm in enumerate(["+", "-"]):
            self.insert_plusminus(cols_right, "col", pm)
        rows_ = Frame(main_left)
        rows_.pack(side="bottom", fill="y", expand=True)
        for i, pm in enumerate(["+", "-"]):
            self.insert_plusminus(rows_, "row", pm)

        self.cp_frame = Frame(main_right)
        self.load_couples(self.cp_frame)
        self.insert_sendcancel(frame=main_bottom)
        self.cp_frame.pack(side="left", fill="both", expand=False)
        self.window.wait_window()
