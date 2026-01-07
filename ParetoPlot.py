import math
import os
import argparse
import matplotlib
from scipy import stats

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt

class ParetoPlot:
    """
    Given a solved PB it plots its Pareto Front, the Nash and Fair solutions
    """

    def __init__(self, solutions, proposals_data, budget, city_bud, citizens, part, real_solution, theta):
        """
        Initialises the variables of the object

        Parameters:
            solutions (list): A list of solutions, each solution should be represented as a list of its selected proposal ids
            proposals_data (list): A list of proposals, each proposal is represented as a lists of its 3 stakeholder utilities
            budget (int): The PB budget
            city_bud (int): The municipal budget of the city (by definition citybud > budget)
            citizens (int): The total population of the city
            part (int): The number of participants the PB process had
            real_solution (list): A list of the proposal ids of the proposals selected in the actual solution of the PB
            theta (float): The minimum contribution percentage parameter theta for the Nash solution
        """
        self.solutions = solutions
        self.data = proposals_data
        self.budget = budget
        self.city_bud = city_bud
        self.citizens = citizens
        self.part = part
        self.theta = theta
        self.sat_points = []
        self.alg_points = []
        self.npar_points = []
        if real_solution:
            self.real_solution = real_solution
        else:
            self.real_solution = None

        # Calculate the 3 stakeholder utilities for each solution
        for s in self.solutions:
            par, gov, npar = self.utils(s)
            self.sat_points.append(par)
            self.alg_points.append(gov)
            self.npar_points.append(npar)

        # We will need the maximum of the three stakeholder utilities for some calculations
        self.maxsat = max(self.sat_points)
        self.maxal = max(self.alg_points)
        self.maxnpar = max(self.npar_points)

    def plot_pareto_and_sols(self, savefile):
        """
        Plots the whole Pareto front, Nash and Fair solutions, and the ideal point.
        The plots have normalised axis running from 0 to 1 where 1 represents the maximum possible
        attainable utility for one of the three stakeholders in the MSPB problem at hand.

        Parameters:
             savefile (str): The name of the file to save the plot
        """

        # Create the figure
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        # We make a copy of sat_points, alg_points, and npar_points so that we don't delete them in the process
        sat_points = []
        alg_points = []
        npar_points = []
        for s in self.sat_points:
            sat_points.append(s)
        for a in self.alg_points:
            alg_points.append(a)
        for a in self.npar_points:
            npar_points.append(a)

        # We plot the ideal point, that with 100% of the maximum possible utility for the 3 stakeholders
        ax.scatter([1], [1], [1], c="green",s=170, marker="*", label="Ideal")
        #If there are points in the Pareto front that have not been plotted

        # We find the axis positions for the Nash and Fair solution (that is their percentage of the maximum
        # possible utility for each of the 3 stakeholders
        nash_sol = [None, None, None]
        fair_sol = [None, None, None]
        nash_pos, nash_sol[0], nash_sol[1], nash_sol[2] = self.nash_sol()
        fair_pos, fair_sol[0], fair_sol[1], fair_sol[2] = self.fair_sol()
        for i in range(len(sat_points)):
            sat_points[i] = sat_points[i]/self.maxsat
            alg_points[i] = alg_points[i]/self.maxal
            npar_points[i] = npar_points[i]/self.maxnpar

        # If we are asked to plot the real solution too, we find its axis positions and print the results
        real_sol_utils = None
        if self.real_solution:
            aux = self.utils(self.real_solution)
            real_sol_utils = [aux[0], aux[1], aux[2]]
            real_sol_utils[0] = real_sol_utils[0] / self.maxsat
            real_sol_utils[1] = real_sol_utils[1] / self.maxal
            real_sol_utils[2] = real_sol_utils[2] / self.maxnpar
            print("The real solution got:"+str(round(real_sol_utils[0]*100,2))+"% of the max votes, "+str(round(real_sol_utils[1]*100,2))+"% of the max govalign, "+str(round(real_sol_utils[2]*100,2))+"% of the max nparlign)")
            ax.scatter(real_sol_utils[0], real_sol_utils[1], real_sol_utils[2], c="black",s=170, marker="X", label="Real solution")

        # We print the results for the Nash and Fair solutions
        print("The Nash solution gets: "+str(round(nash_sol[0]*100,2))+"% of the max votes, "+str(round(nash_sol[1]*100,2))+"% of the max govalign, "+str(round(nash_sol[2]*100,2))+"% of the max nparlign")
        print("The Fair solution gets: "+str(round(fair_sol[0]*100,2))+"% of the max votes, "+str(round(fair_sol[1]*100,2))+"% of the max govalign, "+str(round(fair_sol[2]*100,2))+"% of the max nparlign")

        # In order to plot the Nash and Fair solutions with a different style of point, we remove them from the
        # rest of solutions to be plotted as the Pareto front. If the Nash and Fair solutions are the same, we
        # only plot one point representing both.
        if nash_pos == fair_pos:
            sat_points.pop(nash_pos)
            alg_points.pop(nash_pos)
            npar_points.pop(nash_pos)
        elif nash_pos > fair_pos:
            sat_points.pop(nash_pos)
            alg_points.pop(nash_pos)
            npar_points.pop(nash_pos)
            sat_points.pop(fair_pos)
            alg_points.pop(fair_pos)
            npar_points.pop(fair_pos)
        else:
            sat_points.pop(fair_pos)
            alg_points.pop(fair_pos)
            npar_points.pop(fair_pos)
            sat_points.pop(nash_pos)
            alg_points.pop(nash_pos)
            npar_points.pop(nash_pos)

        # If there are points in the Pareto front after having removed the Nash and Fair solutions we plot them
        if sat_points:
            ax.scatter(sat_points, alg_points, npar_points, c="blue", label="Pareto front", zorder = 1)

        # We plot the Nash and Fair solutions in a different color, if they coincide we only plot one point
        if nash_pos == fair_pos:
            ax.scatter(nash_sol[0], nash_sol[1], nash_sol[2], c="orange",s=170, marker="P", label="Nash and Fair solutions", zorder=2)
        else:
            ax.scatter(nash_sol[0], nash_sol[1], nash_sol[2], c="red",s=170, marker="P", label="Nash solution", zorder=2)
            ax.scatter(fair_sol[0], fair_sol[1], fair_sol[2], c="orange",s=110, marker="D", label="Fair solution", zorder=3)

        # Add the axis labels and the legend
        plt.xlabel("Participant votes")
        plt.ylabel("Government alignment")
        ax.set_zlabel("Non-participant alignment", rotation=90)
        ax.zaxis.labelpad = 10
        plt.legend(loc = "center left")
        fig.tight_layout(pad=0)

        # Save the plot in a new pdf file
        if self.real_solution:
            newsavefile = savefile
            newsavefile = newsavefile.replace(".pdf", "_real.pdf")
            plt.savefig(newsavefile)
        else:
            plt.savefig(savefile)
        plt.close()

        # We now plot and save the pdfs for the three different 2d projections
        self.two_stakeholder_plot(sat_points, alg_points, npar_points, nash_sol, fair_sol, savefile, real_sol_utils)

    def contributions(self):
        """
            Calculates the contribution for each stakeholder, according to the formula in our paper
            :return: the 3 contributions as floats
        """

        # citizen budget share, the amount of budget contributed by each citizen (through taxes)
        cbs = self.city_bud / self.citizens

        # Contributions for participants, government, and non-participants
        cpc = min(cbs * self.part, self.city_bud - self.theta * self.city_bud)
        cgov = self.budget
        cnp = min(cbs * (self.citizens - self.part), self.city_bud - self.theta * self.city_bud)

        """
            We have to ensure the contributions are above the minimum percentage threshold (the theta parameter)
            Therefore we have to check 6 cases (2 cases where 2 of them are below theta, 3 cases where 1 of them
            is below theta, and 1 case where none of them are below theta. The other two cases cannot happen,
            this is because the contribution for participants and non-participants are complementary, therefore 
            it is impossible that both of them are below theta, and therefore it is also impossible that all 
            contributions are below theta.
        """

        # If the participants contribution is less than the minimum
        if cpc < self.theta * self.city_bud:

            # Case 1: Participant and government contribution are less than the minimum (theta)
            # We assign theta to both of them and adjust the non-participant contribution
            if cgov < self.theta * self.city_bud:
                cpc = self.theta
                cgov = self.theta
                cnp = 1 - 2 * self.theta

            # Case 2: If only the participants contribution is less than the minimum, we rescale the
            # rest of contributions
            else:
                cpc = self.theta * self.city_bud
                tot = cpc + cnp + cgov
                cnp = cnp / tot
                cgov = cgov / tot

        # If the government contribution is less than the minimum
        elif cgov < self.theta * self.city_bud:

            # Case 3: Government and non-participant contribution are below theta. We assign the
            # minimum and rescale the participants contribution
            if cnp < self.theta * self.city_bud:
                cgov = self.theta
                cnp = self.theta
                cpc = 1 - 2 * self.theta

            # Case 4: Only the government contribution is below theta.
            else:
                cgov = self.theta * self.city_bud
                tot = cpc + cgov + cnp
                cpc = cpc / tot
                cnp = cnp / tot

        # Case 5: Non-participants contribution is the only one less than the minimum
        # We assign the minimum to cnp and rescale the other contributions
        elif cnp < self.theta * self.city_bud:
            cnp = self.theta * self.city_bud
            tot = cpc + cgov + cnp
            cpc = cgov / tot
            cgov = cgov / tot

        # Case 6: If no contribution is below the minimum, we just normalise them
        else:
            tot = cpc + cgov + cnp
            cpc = cpc / tot
            cgov = cgov / tot
            cnp = cnp / tot

        return cpc, cgov, cnp

    def nash_sol(self):
        """
        Find the Nash solution of those in the Pareto front.

        Returns:
            nas_sol_pos (int): The solution id of the Nash solution
            nas_sat (float): The citizen satisfaction of the Nash solution.
            nash_al (float): The government alignment of the Nash solution.
            nash_npar (float): The non-participant alignment of the Nash solution.
        """

        # Initialise the Nash product and the solution id of the Nash solution
        nash_product = 0
        nash_sol_pos = None

        #Calculate the contributions
        cpc, cgov, cnp = self.contributions()

        # Print the contributions
        print("The contributions are: "+str(cpc)+" (participants),"+str(cgov)+" (government),"+str(cnp)+" (non-participants)")

        # For each solution in the Pareto front, calculate the nash product
        for i in range(len(self.sat_points)):
            prod = ((self.sat_points[i]/self.maxsat) ** cpc) * ((self.alg_points[i]/self.maxal) ** cgov) * ((self.npar_points[i]/self.maxnpar) ** cnp)

            # If the product is greater than that of the previous maximum, update the maximum and the Nash sol id
            if prod > nash_product:
                nash_product = prod
                nash_sol_pos = i

        # Calculate the plot points for the Nash solution (percentages over the maximum for each utility)
        nash_sat = self.sat_points[nash_sol_pos] / self.maxsat
        nash_al = self.alg_points[nash_sol_pos] / self.maxal
        nash_npar = self.npar_points[nash_sol_pos] / self.maxnpar

        return nash_sol_pos, nash_sat, nash_al , nash_npar

    def fair_sol(self):
        """
        Calculate the Fair solution of those in the Pareto front.

        Returns:
            fair_sol_pos (int): The solution id of the Fair solution
            fair_sat (float): The citizen satisfaction of the Fair solution.
            fair_gov (float): The government alignment of the Fair solution.
            fair_npar (float): The non-participant alignment of the Fair solution.
        """

        # Find the 3-stakeholder Nash solution

        nash_sol_id = self.nash_sol()[0]

        #Calculate the contributions of each stakeholder needed for the worth of 2-stakeholder coalitions
        cpc, cgov, cnp = self.contributions()

        # Find the 2-stakeholder Nash solution for participants and government
        nash_product = 0
        nash_pos_par_gov = None
        if cpc + cgov > 0.5:
            for i in range(len(self.sat_points)):
                prod = (self.sat_points[i]/self.maxsat) * (self.alg_points[i]/self.maxal)
                if prod > nash_product:
                    nash_product = prod
                    nash_pos_par_gov = i

        # Find the 2-stakeholder Nash solution for participants and non-participants
        nash_product = 0
        nash_pos_par_npar = None
        if cpc + cnp > 0.5:
            for i in range(len(self.sat_points)):
                prod = (self.sat_points[i]/self.maxsat) * (self.npar_points[i]/self.maxnpar)
                if prod > nash_product:
                    nash_product = prod
                    nash_pos_par_npar = i

        # Find the 2-stakeholder Nash solution for government and non-participants
        nash_product = 0
        nash_pos_gov_npar = None
        if cgov + cnp > 0.5:
            for i in range(len(self.alg_points)):
                prod = (self.alg_points[i]/self.maxal) * (self.npar_points[i]/self.maxnpar)
                if prod > nash_product:
                    nash_product = prod
                    nash_pos_gov_npar = i

        # Initialise the fair solution id and the minimum complaint to find the nucleolus
        fair_sol_pos = None
        min_complaint = [math.inf]*7

        # For each solution we calculate the complaint vector for the stakeholder coalitions
        for i in range(len(self.sat_points)):

            # We initialise the complaint vector with a 0 (the complaint of the empty coalition)
            complaint = [0]

            #We calculate the complaint with regards to the singleton coalitions (their worth is always 0)
            # Complaint for the participant singleton
            c = -self.sat_points[i]/self.maxsat
            complaint.append(c)

            # Complaint for the government singleton
            c = -self.alg_points[i]/self.maxal
            complaint.append(c)

            # Complaint for non-participant singleton
            c = -self.npar_points[i]/self.maxnpar
            complaint.append(c)

            # Complaint with regard to the participant and government Nash solution
            if nash_pos_par_gov:
                c = (self.sat_points[nash_pos_par_gov]-self.sat_points[i])/self.maxsat
                c += (self.alg_points[nash_pos_par_gov]-self.alg_points[i])/self.maxal
            else:
                c = -self.sat_points[i] / self.maxsat
                c -= self.alg_points[i]/self.maxal
            complaint.append(c)

            # Complaint with regard to the participant and non-participant Nash solution
            if nash_pos_par_npar:
                c = (self.sat_points[nash_pos_par_npar] - self.sat_points[i])/self.maxsat
                c += (self.npar_points[nash_pos_par_npar] - self.npar_points[i])/self.maxnpar
            else:
                c = -self.sat_points[i] / self.maxsat
                c -= self.npar_points[i] / self.maxnpar
            complaint.append(c)

            # Disparity with regard to the government and non-participant Nash solution
            if nash_pos_gov_npar:
                c = (self.alg_points[nash_pos_gov_npar] - self.alg_points[i])/self.maxal
                c += (self.npar_points[nash_pos_gov_npar] - self.npar_points[i])/self.maxnpar
            else:
                c = -self.alg_points[i] / self.maxal
                c -= self.npar_points[i] / self.maxnpar
            complaint.append(c)

            #Camplaint for the 3-stakeholder coalition
            c = (self.sat_points[nash_sol_id]-self.sat_points[i])/self.maxsat
            c += (self.alg_points[nash_sol_id]-self.alg_points[i])/self.maxal
            c += (self.npar_points[nash_sol_id]-self.npar_points[i])/self.maxnpar
            complaint.append(c)

            #We order the complaint vector from maximum complaint to minimum
            complaint.sort(reverse=True)

            # The current solution is the new Fair solution if its complaint vector is lexicographically smaller than
            # that of our previous Fair solution
            if complaint < min_complaint:
                min_complaint = complaint.copy()
                fair_sol_pos = i

        # Calculate the plot points for the fair solution (percentages over the maximum for each utility)
        fair_sat = self.sat_points[fair_sol_pos] / self.maxsat
        fair_gov = self.alg_points[fair_sol_pos] / self.maxal
        fair_npar = self.npar_points[fair_sol_pos] / self.maxnpar

        return fair_sol_pos, fair_sat, fair_gov, fair_npar

    def utils(self, proposals):
        """
        Calculates the citizen satisfaction, government alignment, and non-participant alignment of a set of proposals.
        Normally representing a solution.

        Parameters:
            proposals (list): list of proposals

        Returns:
            sat (float): citizen satisfaction of the proposals
            goval (float): government alignment of the proposals
            nparal (float): non-participant alignment of the proposals
        """
        sat = 0
        goval = 0
        nparal = 0
        for p in proposals:
            sat += self.data[p-1][0]
            goval += self.data[p-1][1]
            nparal += self.data[p-1][2]
        return sat, goval, nparal

    def two_stakeholder_plot(self, satpoints, algpoints, nparpoints, nashsol, fairsol, savefile, realsol = None):
        """
        Plots the 2D projections of the Pareto Front and saves then in PDF files.

        Parameters:
            satpoints (list): list of the citizen satisfaction utilities for each proposal id
            algpoints (list): list of the government alignment utilities for each proposal id
            nparpoints (list): list of the non-participant alignment utilities for each proposal id
            nashsol (list): the 3 stakeholder utilities of the Nash solution
            fairsol (list): the 3 stakeholder utilities of the Fair solution
            savefile (string): the name of the file to save the PDF
            realsol (list): the 3 stakeholder utilities of the real solution
        """
        newsavefile = savefile.replace(".pdf", "_pargov.pdf")
        self.par_gov_plot(satpoints, algpoints, nashsol, fairsol, newsavefile, realsol)
        newsavefile = savefile.replace(".pdf", "_parnpar.pdf")
        self.par_npar_plot(satpoints, nparpoints, nashsol, fairsol, newsavefile, realsol)
        newsavefile = savefile.replace(".pdf", "_govnpar.pdf")
        self.gov_npar_plot(algpoints, nparpoints, nashsol, fairsol, newsavefile, realsol)

    def par_gov_plot(self, satpoints, algpoints, nash_sol, fair_sol, savefile, real_sol = None):
        """
        Plots the 2D projection of the Pareto Front over the plane formed by the citizen satisfaction and government
        alignment axis.

        Parameters:
            satpoints (list): list of the citizen satisfaction utilities for each proposal id
            algpoints (list): list of the government alignment utilities for each proposal id
            nash_sol (list): the 3 stakeholder utilities of the Nash solution
            fair_sol (list): the 3 stakeholder utilities of the Fair solution
            savefile (string): the name of the file to save the PDF
            real_sol (list): the 3 stakeholder utilities of the real solution
        """

        # Create the figure
        fig = plt.figure()
        ax = fig.add_subplot()

        # If there are points in the Pareto front other than our solutions, plot the points of the Pareto front
        # and print stats for our experiment discussion.
        if satpoints:
            ax.scatter(satpoints, algpoints, c="blue", label="Pareto front", zorder=1)
            coeff = stats.pearsonr(satpoints, algpoints)
            print("Par Gov Pearson coeff:" + str(coeff.statistic) + " (p-value: " + str(coeff.pvalue) + ")")
            coeff = stats.spearmanr(satpoints, algpoints)
            print("Par Gov Spearman coeff:"+str(coeff.statistic)+" (p-value: "+str(coeff.pvalue)+")")
            coeff = stats.kendalltau(satpoints, algpoints)
            print("Par Gov Kendall coeff:" + str(coeff.statistic) + " (p-value: " + str(coeff.pvalue) + ")")

        # If required, plot the point of the real solution's utilities
        if real_sol:
            ax.scatter(real_sol[0], real_sol[1], c="black", s=170, marker="X", label="Real solution", zorder=2)

        # Plot the Nash and Fair solutions (if they are the same plot only one point)
        if nash_sol == fair_sol:
            ax.scatter(nash_sol[0], nash_sol[1], c="orange", s=170, marker="P", label="Nash and Fair solution", zorder=3)
        else:
            ax.scatter(nash_sol[0], nash_sol[1], c="red", s=170, marker="P", label="Nash solution", zorder=3)
            ax.scatter(fair_sol[0], fair_sol[1], c="orange", s=110, marker="D", label="Fair solution", zorder=4)

        # Add the axis labels and the legend
        plt.xlabel("Participant votes")
        plt.ylabel("Government alignment")
        plt.legend()

        # Save the figure in a PDF file
        plt.savefig(savefile,bbox_inches='tight')

    def par_npar_plot(self, satpoints, nparpoints, nash_sol, fair_sol, savefile, realsol = None):
        """
        Plots the 2D projection of the Pareto Front over the plane formed by the citizen satisfaction and
        non-participant alignment axis.

        Parameters:
            satpoints (list): list of the citizen satisfaction utilities for each proposal id
            nparpoints (list): list of the non-participant alignment utilities for each proposal id
            nash_sol (list): the 3 stakeholder utilities of the Nash solution
            fair_sol (list): the 3 stakeholder utilities of the Fair solution
            savefile (string): the name of the file to save the PDF
            realsol (list): the 3 stakeholder utilities of the real solution
        """

        # Create the figure
        fig = plt.figure()
        ax = fig.add_subplot()

        # If there are points in the Pareto front other than our solutions, plot the points of the Pareto front
        # and print stats for our experiment discussion.
        if satpoints:
            ax.scatter(satpoints, nparpoints, c="blue", label="Pareto front", zorder=1)
            coeff = stats.pearsonr(satpoints, nparpoints)
            print("Par NPar Pearson coeff:" + str(coeff.statistic) + " (p-value: " + str(coeff.pvalue) + ")")
            coeff = stats.spearmanr(satpoints, nparpoints)
            print("Par NPar Spearman coeff:" + str(coeff.statistic) + " (p-value: " + str(coeff.pvalue) + ")")
            coeff = stats.kendalltau(satpoints, nparpoints)
            print("Par NPar Kendall coeff:" + str(coeff.statistic) + " (p-value: " + str(coeff.pvalue) + ")")

        # If required, plot the point of the real solution's utilities
        if realsol:
            ax.scatter(realsol[0], realsol[2], c="black", s=170, marker="X", label="Real solution", zorder=2)

        # Plot the Nash and Fair solutions (if they are the same plot only one point)
        if nash_sol == fair_sol:
            ax.scatter(nash_sol[0], nash_sol[2], c="orange", s=170, marker="P", label="Nash and Fair solution", zorder=3)
        else:
            ax.scatter(nash_sol[0], nash_sol[2], c="red", s=170, marker="P", label="Nash solution", zorder=3)
            ax.scatter(fair_sol[0], fair_sol[2], c="orange", s=110, marker="D", label="Fair solution", zorder=4)

        # Add the axis labels and the legend
        plt.xlabel("Participant votes")
        plt.ylabel("Non-participant alignment")
        plt.legend()

        # Save the figure in a PDF file
        plt.savefig(savefile,bbox_inches='tight')

    def gov_npar_plot(self, algpoints, nparpoints, nash_sol, fair_sol, savefile, realsol = None):
        """
        Plots the 2D projection of the Pareto Front over the plane formed by the government alignment and
         non-participant alignment axis.

        Parameters:
            algpoints (list): list of the government alignment utilities for each proposal id
            nparpoints (list): list of the non-participant alignment utilities for each proposal id
            nash_sol (list): the 3 stakeholder utilities of the Nash solution
            fair_sol (list): the 3 stakeholder utilities of the Fair solution
            savefile (string): the name of the file to save the PDF
            realsol (list): the 3 stakeholder utilities of the real solution
        """

        # Create the figure
        fig = plt.figure()
        ax = fig.add_subplot()

        # If there are points in the Pareto front other than our solutions, plot the points of the Pareto front
        # and print stats for our experiment discussion.
        if algpoints:
            ax.scatter(algpoints, nparpoints, c="blue", label="Pareto front", zorder=1)
            coeff = stats.pearsonr(algpoints, nparpoints)
            print("Gov NPar Pearson coeff:" + str(coeff.statistic) + " (p-value: " + str(coeff.pvalue) + ")")
            coeff = stats.spearmanr(algpoints, nparpoints)
            print("Gov NPar Spearman coeff:" + str(coeff.statistic) + " (p-value: " + str(coeff.pvalue) + ")")
            coeff = stats.kendalltau(algpoints, nparpoints)
            print("Gov NPar Kendall coeff:" + str(coeff.statistic) + " (p-value: " + str(coeff.pvalue) + ")")

        # If required, plot the point of the real solution's utilities
        if realsol:
            ax.scatter(realsol[1], realsol[2], c="black", s=170, marker="X", label="Real solution", zorder=2)

        # Plot the Nash and Fair solutions (if they are the same plot only one point)
        if nash_sol == fair_sol:
            ax.scatter(nash_sol[1], nash_sol[2], c="orange", s=170, marker="P", label="Nash and Fair solution", zorder=3)
        else:
            ax.scatter(nash_sol[1], nash_sol[2], c="red", s=170, marker="P", label="Nash solution", zorder=3)
            ax.scatter(fair_sol[1], fair_sol[2], c="orange", s=110, marker="D", label="Fair solution", zorder=4)

        # Add the axis labels and legend to the plot
        plt.xlabel("Government alignment")
        plt.ylabel("Non-participant alignment")
        plt.legend()

        # Save the figure in a PDF file
        plt.savefig(savefile,bbox_inches='tight')

def main():

    # We process the --pb --theta and --real_sol arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--pb", required=True, type=str)
    parser.add_argument("--theta", required=False, type=float, default=0.2)
    parser.add_argument("--real_sol", required=False, action = "store_true")
    args = parser.parse_args()
    pb = args.pb
    theta = args.theta
    plot_real_sol = args.real_sol

    # We check the arguments have been provided correctly
    if pb in ["Barcelona", "Warszawa23", "Warszawa24"] and 0 <= theta <= 1:

        # Read the file with the basic PB data, utilities for each project, and the actual solution of the PB
        juliasource = open(os.getcwd() + "/JuliaDataFiles/" + pb + ".txt", "r")
        lines = juliasource.readlines()
        juliasource.close()

        # Store 3 stakeholder utilities for each project in the PB
        citsat = lines[0].split("=")[1].replace(" ", "").replace("\n", "").replace("[", "").replace("]", "").split(",")
        govalignment = lines[1].split("=")[1].replace(" ", "").replace("\n", "").replace("[", "").replace("]", "").split(",")
        nparalignment = lines[2].split("=")[1].replace(" ", "").replace("\n", "").replace("[", "").replace("]", "").split(",")
        proposalsdata = []
        for i in range(len(citsat)):
            proposalsdata.append([int(citsat[i]), int(govalignment[i]), int(nparalignment[i])])

        # Store the basic PB information
        budget = int(lines[4].split("=")[1].replace("\n", "").replace(" ", ""))
        citybudget = int(lines[5].split("=")[1].replace("\n", "").replace(" ", ""))
        population = int(lines[6].split("=")[1].replace("\n", "").replace(" ", ""))
        participants = int(lines[7].split("=")[1].replace("\n", "").replace(" ", ""))

        # Store the projects that were selected in the actual PB process if required
        if plot_real_sol:
            real_solution = []
            selected = lines[8].split("=")[1].replace(" ", "").replace("\n", "").replace("[", "").replace("]", "").split(",")
            for i in selected:
                real_solution.append(int(i))
        else:
            real_solution = None

        # Now we read the file containing the data of the solutions in the Pareto front of the MSPB problem
        datasource = open(os.getcwd()+"/results/data/"+"PB_"+pb+".txt", "r")
        lines = datasource.readlines()
        datasource.close()

        # The Julia packages we use to calculate the Pareto front sometimes include repeated solutions
        # (which can be problematic down the line), so we check this to ensure everything is fine
        solproposals = []
        isproposal = False
        count = 0
        linenum = 0
        repeatedsols = []
        for l in lines[1:]:
            if "list of selected proposals" in l:
                isproposal = True
            if not l:
                isproposal = False

            if isproposal and not "list of selected proposals" in l:
                linenum += 1
                l = l.split(",")
                sol_props = []
                for p in l[1:]:
                    sol_props.append(int(p))
                sol_props.sort()
                if sol_props in solproposals:
                    count += 1
                    repeatedsols.append(linenum-1)
                else:
                    solproposals.append(sol_props)

        # We print information on repeated solutions for debugging purposes
        print("Number of repeated solutions:"+str(count))
        print("Repeated solutions:"+str(repeatedsols))

        # We instantiate a ParetoPlot object with all the information
        pareto = ParetoPlot(solproposals, proposalsdata, budget, citybudget, population, participants, real_solution, theta)

        # We plot the points and save them in the results/plots directory
        filename = os.getcwd()+"/results/plots/"+pb+"/"+pb+".pdf"
        pareto.plot_pareto_and_sols(filename)

    # If arguments were provided incorrectly we print error messages
    elif 0 <= theta <= 1./3.:
        print("Invalid PB process. Valid PBs are: Barcelona, Warszawa23, and Warszawa24")
    else:
        print("Invalid theta value. Theta has to be a float between 0 and 1.")

if __name__ == "__main__":
    main()