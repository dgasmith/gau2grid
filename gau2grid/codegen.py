"""
Contains a simple codegen helper
"""


class CodeGen(object):
    def __init__(self, indent="    ", cgen=False):
        self.indent_lvl = 0
        self.indent_tab = indent
        self.data = []
        self.cgen = cgen

    def indent(self, lvl=1):
        """
        Indents the code one or more levels
        """

        self.indent_lvl += lvl

    def dedent(self, lvl=1):
        """
        Dedents the code one or more levels
        """

        self.indent_lvl -= lvl
        if self.indent_lvl < 0:
            last_lines = "\n".join(self.data[-4:])
            raise ValueError("Indent level is negative! Last lines:\n\n%s" % last_lines)

    def write(self, line, endl=None):
        """
        Write a line with the current indent
        """
        shift = self.indent_lvl * self.indent_tab
        if self.cgen and (endl is None) and ("//" not in line) and ("#" not in line):
            endl = ";"
        if endl is None:
            endl = ""

        self.data.append(shift + line + endl)

    def blankline(self):
        """
        Inserts a blankline
        """
        self.data.append("")

    def repr(self, filename=None, combine="\n", clang_format=False):
        """
        Combined the data into a single string
        """
        tmp = combine.join(self.data)
        if clang_format:
            if self.cgen is False:
                raise KeyError("clang_format is only valid for c generation.")
            try:
                tmp = run_clang_format(tmp)
            except:
                print("Could not run clang-format, skipping")

        if filename is None:
            return tmp
        else:
            with open(filename, "w") as outfile:
                outfile.write(tmp)

    def start_c_block(self, line=None):
        """
        Opens a C block with open brackets and indention
        """

        if self.cgen is False:
            raise KeyError("Start c block only valid for c generation.")

        if line is None:
            self.write("{", endl="")
        else:
            self.write(line + " {", endl="")

        self.indent()

    def close_c_block(self):
        """
        Ends a c block with a dedent and close line
        """
        if self.cgen is False:
            raise KeyError("Start c block only valid for c generation.")

        self.dedent()
        self.write("}", endl="")


def run_clang_format(text):
    import subprocess as sp

    # Gotta escape those characters
    text = text.replace('"', '\\"')

    cmd = 'echo "' + text + '"'
    cmd += " | clang-format"
    output = sp.check_output(cmd, shell=True)
    return output.decode()
