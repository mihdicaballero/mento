from IPython.display import Markdown, display
import pandas as pd

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
        if self.mid_value > FU:
            color = self.green
        elif self.mid_value <= FU <= self.max_value:
            color = self.yellow
        else:
            color = self.red
        
        return f"$\\color{{{color}}}{{\\text{{FU}}={round(FU,2)}}}$"
    def FU_value(self, FU):
        # Determine color based on FU value
        if self.mid_value > FU:
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
            if self.mid_value > FU:
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
    

class SectionSummary:
    def __init__(self, sections):
        """
        Initializes the SectionSummary with a list of Section objects.
        
        :param sections: List of Section objects
        """
        self.capacity_columns = ["Viga", "b", "h", "As.inf", "As.sup", "Av", 
                        "As.inf.real", "As.sup.real", "Av.real", "ØMn", "ØVn"]
        self.check_columns = ["Viga", "b", "h", "Vu", "Mu", "As.inf.real", "As.inf.nec", 
                         "As.sup.real", "As.sup.nec", "Av.real", "Av.nec", "MRd.inf","MRd.sup", 
                         "VRd", "FUb.inf", "FUb.sup", "FUv"]
        self.sections = sections  # Store the list of sections

    def capacity(self):
        """
        Creates a DataFrame from the list of Section objects and adds the units row.
        
        :return: A pandas DataFrame with the data from Section objects and units as the first row.
        """
        # Convert each section object to a row dictionary
        data = [section.to_capacity_dict() for section in self.sections]
        
        # Add units row
        units_row = {
            "Viga": "",
            "b": "cm",
            "h": "cm",
            "As.inf": "",
            "As.sup": "",
            "Av": "",
            "As.inf.real": "cm²",
            "As.sup.real": "cm²",
            "Av.real": "cm²/m",
            "MRd": "kNm",
            "VRd": "kN",
        }
        
        # Combine units row with section data
        data.insert(0, units_row)
        
        # Create a DataFrame from the data
        df = pd.DataFrame(data, columns=self.capacity_columns)
        
        return df
    def check(self):
        """
        Creates a DataFrame for checking shear and bending for each section.
        Outputs necessary and real reinforcement, utilization factors for bending and shear.
        
        :return: A pandas DataFrame summarizing the check of each section.
        """

        # Convert each section object to a row dictionary
        data = [section.to_check_dict() for section in self.sections]

        # Create a DataFrame from the data
        df = pd.DataFrame(data, columns=self.check_columns)

        # Add units row
        units_row = {
            "Viga": "",
            "b": "cm",
            "h": "cm",
            "Vu": "kN",
            "Mu": "kNm",
            "As.inf.real": "cm²",
            "As.inf.nec": "cm²",
            "As.sup.real": "cm²",
            "As.sup.nec": "cm²",
            "Av.real": "cm²/m",
            "Av.nec": "cm²/m",
            "MRd.inf": "kNm",
            "MRd.sup": "kNm",
            "VRd": "kN",
            "FUb.inf": "", 
            "FUb.sup": "",
            "FUv": "" 
        }
        
        # Append the units row to the top of the DataFrame
        df = pd.concat([pd.DataFrame([units_row]), df], ignore_index=True)
        
        # Add style to FU columns
        self.fu_columns = ["FUb.inf", "FUb.sup", "FUv"]  # Columns to be styled based on FU value
        formatter = Formatter()
        # Apply color formatting to the specified FU columns
        return formatter.color_FU_df(df, self.fu_columns)

def main():
    #  Examples to run in a Jupyter Notebook
    formatter = Formatter()
    display(Markdown(formatter.FU(0.85)))
    display(Markdown(formatter.FU_value(0.85)))
    display(Markdown(formatter.is_lower(0.85,1)))
    display(Markdown(formatter.is_greater(0.85,1)))

if __name__ == "__main__":
    main()