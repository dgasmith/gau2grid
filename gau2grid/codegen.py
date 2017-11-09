"""
Contains a simple codegen helper
"""


class CodeGen(object):
    def __init__(self, indent="    "):
        self.indent_lvl = 0
        self.indent_tab = indent
        self.data = []

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
            raise ValueError("Indent level is negative!")

    def write(self, line):
        """
        Write a line with the current indent
        """
        shift = self.indent_lvl * self.indent_tab
        self.data.append(shift + line)

    def blankline(self):
        """
        Inserts a blankline
        """
        self.data.append("")

    def repr(self, combine="\n"):
        """
        Combined the data into a single string
        """
        return combine.join(self.data)

    def start_c_block(self, line=None):
        """
        Opens a C block with open brackets and indention
        """

        if line is None:
            self.write("{")
        else:
            self.write(line + " {")

        self.indent()

    def stop_c_block():
        """
        Ends a c block with a dedent and close line
        """

        self.dedent()
        self.write("}")



