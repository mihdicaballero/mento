from IPython.display import Markdown, display

class Formatter:
    """
    Class to stly ethe results in a Markdown display and a DataFrame output.
    """
    def __init__(self):
        # Define colors in the constructor
        self.green = '#439b00'
        self.red = '#d43e36'
        self.yellow = '#efc200'
        self.mid_value=0.95
        self.max_value=1
    
    def FU(self, FU):
        # Determine color based on FU value
        if FU < self.mid_value:
            color = self.green
        elif self.mid_value <= FU <= self.max_value:
            color = self.yellow
        else:
            color = self.red
        
        return f"$\\color{{{color}}}{{\\text{{FU}}={round(FU,2)}}}$"
    def FU_value(self, FU):
        # Determine color based on FU value
        if FU < self.mid_value:
            color = self.green
        elif self.mid_value <= FU <= 1:
            color = self.yellow
        else:
            color = self.red
        
        return f"$\\color{{{color}}}{{{round(FU,2)}}}$"
    
    def is_lower(self, value1, value2):
        # Compare two values and return the appropriate formatted output
        if value1 < value2:
            return r"$\color{" + self.green + r"}{\, \checkmark}$"  # Green checkmark
        else:
            return r"$\color{" + self.red + r"}{\, \times}$"  # Red cross
    
    def is_greater(self, value1, value2):
        # Compare two values and return the appropriate formatted output
        if value1 > value2:
            return r"$\color{" + self.green + r"}{\, \checkmark}$"  # Green checkmark
        else:
            return r"$\color{" + self.red + r"}{\, \times}$"  # Red cross
    # Formatting functions outside the class
    def FU_value_df(self,FU):
        if isinstance(FU, (int, float)):
            if FU < self.mid_value:
                return f"color: {self.green}"
            elif self.mid_value <= FU <= self.max_value:
                return f"color: {self.yellow}"
            else:
                return f"color: {self.red}"
        else: 
            return ""

    def apply_FU_style(self,value):
        formatted_value = round(value, 2) if isinstance(value, (int, float)) else value
        return self.FU_value_df(formatted_value)

    def color_FU_df(self,df, fu_columns):
        """
        Apply color styling to specified FU-related columns in the DataFrame.
        
        :param df: DataFrame with FU values to color.
        :param fu_columns: List of column names to apply the FU styling to.
        :return: A styled DataFrame with colored FU values.
        """
        return df.style.applymap(self.apply_FU_style, subset=fu_columns).format(precision=2)
    

def main():
    #  Examples to run in a Jupyter Notebook
    formatter = Formatter()
    display(Markdown(formatter.FU(0.85)))
    display(Markdown(formatter.FU_value(0.85)))
    display(Markdown(formatter.is_lower(0.85,1)))
    display(Markdown(formatter.is_greater(0.85,1)))

if __name__ == "__main__":
    main()