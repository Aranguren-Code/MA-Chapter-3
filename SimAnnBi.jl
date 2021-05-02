# Loading the necessary packages
#-------------------------------------------------------------------------------
using CSV
using JuMP
using CPLEX
using DataFrames
using DelimitedFiles

# Load Model ---------------------------------------------
first = 0 # 1 if first run, 0 if you'd like to reset model, 2 if you'd like to continue with previous sol

#---------------------------------------------------------
if first == 0
	print("Loading Model \r\n")
	include("getdata.jl");
	include("Load_Bi.jl");
elseif first == 0  
	print("Reseting Model \r\n")
	include("reset_SA.jl"); #check if correct
end
#-----------------------------------------------

#swaps = open("swaps.csv", "w")
Solutions = open("SA Solutions.csv", "w")
#SA_bases = open("SA_bases.csv", "w")
#write(swaps, "m, s1, id1, s0, id2, \r\n")
write(Solutions,"step,temp, next_energy, delta_e, Best, GHG\r\n")
steps = 20
count_limit = 5
#---------------------------------------------------

start = time() 

let
	current_energy = (maximum(ShortCost) * sum(Demand)) #third party supplier cost 
	Best = current_energy #initial solution
	next_energy = 0 #initialize variable
	
	base = zeros(Int64, D, 2) #ID and state
	
	for i in 1:D
		base[i,1] = i 
		base[i,2] = 0
	end

	z_base = base
	
for g in 1:20 #emissions steps
	
	# current_energy = (maximum(ShortCost) * sum(Demand)) #third party supplier cost 
	# Best = current_energy #initial solution
	# 
	count = 0

	#-----------------------------------------------------------------------
	i_temp = 9000 #TUNING REQUIRED
	f_temp = 4000
	step_t = 250
	ite = Iterators.reverse(f_temp:step_t:i_temp)
	factor = 0.10
	#Change Emission Limit
	JuMP.fix(Epoint, e_step[g])
	print("Emission step $(e_step[g]) \r\n")
	
	#-----------------------------------------------------------------------
	for m in ite 

	print("Running $g Emissions Step at $m Temp \r\n")

		sum_base = sum(base[:,2])
		no_depots =  Int(ceil(rand()*(D*factor))) #to be changed
		delta_d = no_depots - sum_base 

		previous_base = base #Before we make any move

		if delta_d > 0 #Addition
			print("Adding Depots\r\n")
			op_id = 1
			for j in 1:delta_d
				r1 = Int(ceil(rand()*(D-sum_base))) # IF ALL DEPOTS SELECTED, THIS MAY RESULT IN ERROR
				print("$(r1) \r\n")
				base[r1,2] = 1 #change base
				JuMP.fix(D_Model[base[r1,1]], 1) #change model
			end
		elseif delta_d < 0 #Remove
			print("Removing Depots")
			op_id = 0
			delta_d = abs(delta_d)
			for i in 1:delta_d
				r1 = Int(ceil(rand()*sum_base)) + (D-sum_base)
				base[r1,2] = 0
				JuMP.fix(D_Model[base[r1,1]], 0) #check if working
			end
		else #delta_d = 0 and we SWAP
			op_id = 2
			print("Swapping Depots")
			for i in 1:m
				s1 = Int(ceil(rand()*sum_base)) + (D-sum_base)
				s0 = Int(ceil(rand()*(D-sum_base)))
				b1 = base[s1,2]
				b0 = base[s0,2]
				base[s1,2] = b0
				base[s0,2] = b1
				JuMP.fix(D_Model[base[s1,1]], b0)
				JuMP.fix(D_Model[base[s0,1]], b1)
				#write(swaps, "$m, $s1 ,$(base[s1,1]), $s0, $(base[s0,1]) \r\n")
			end	
		end

		#write(SA_bases, "Temp:, $m, Delta_D,$(delta_d),$(op_id),\r\n")
		#writedlm(SA_bases, base, ',')

		#print("initial")
		#print(m)
		#print(base)
#-------------------------------------------------------------------------------------

		optimize!(HS_P3)

#-----------------------------------------------------------------------------------
# Energies and acceptance criteria for simulated annealing
		next_energy = JuMP.objective_value(HS_P3)
		delta_e = current_energy - next_energy
		#write(SA_bases, "Next Energy,$(next_energy),\r\n")
		#write(Solutions,"step,temp, next_energy, delta_e, Best\r\n")
		write(Solutions, "$g,$m, $next_energy, $delta_e, $Best, $(e_step[g]) \r\n")
#-----------------------------------------------------------------------------------------#
#							 Accepting criteria 
#-----------------------------------------------------------------------------------------#
		if delta_e > 0 #accept
			flag = 1
			base = sortslices(base, dims = 1 , by = x -> x[2]) #sort if energy improved
			current_energy = next_energy
#--------------------------------------------------------
			base = z_base #zero out base
			y_avg = 0
			flow_y = zeros(Int64,1)
			for j in 1:D
				for k in 1:C
					for s in 1:S
						if JuMP.value(FLOWY[j,k,s]) > 0
							y_avg += JuMP.value(FLOWY[j,k,s])*prob[s]
						end
					end
					
				end
				if y_avg > 0 
					print(j)
					print(",")
					print(y_avg)
					print("\r\n")
					append!(flow_y,j)
				end
				y_avg = 0
			end
			
			for i in 1:D
				JuMP.fix(D_Model[i],0)
			end

			for i in 2:length(flow_y)
				JuMP.fix(D_Model[flow_y[i]],1)
				base[flow_y[i],2] = 1
			end
				
#--------------------------------------------------------
		elseif delta_e <= 0 
			r_accept = rand()
			if (exp((-D*m)/delta_e)-1) > r_accept #check if met criteria met
					flag = 1 
					current_energy = next_energy
#--------------------------------------------------------
				base = z_base #zero out base
				y_avg = 0
				flow_y = zeros(Int64,1)
				for j in 1:D
					for k in 1:C
						for s in 1:S
							if JuMP.value(FLOWY[j,k,s]) > 0
								y_avg += JuMP.value(FLOWY[j,k,s])*prob[s]
							end
						end
						
					end
					if y_avg > 0 
						print(j)
						print(",")
						print(y_avg)
						print("\r\n")
						append!(flow_y,j)
					end
					y_avg = 0
				end
				
				for i in 1:D
					JuMP.fix(D_Model[i],0)
				end

				for i in 2:length(flow_y)
					JuMP.fix(D_Model[flow_y[i]],1)
					base[flow_y[i],2] = 1
				end
#--------------------------------------------------------
			else 
				base = previous_base

				for i in 1:D
					JuMP.fix(D_Model[i],0)
				end

				for i in 1:D
					if base[i,2] == 1
						JuMP,fix(D_Model[base[i,1]],1)
					end
				end

				flag = 0
				count += 1
				print("Count: $count \r\n")
				y_avg = 0.0
			end
		end

		if next_energy < Best
			count = 0 #count resets
			print("Count Reset due to improved solution \r\n")
			Best = next_energy
			Best_base = base 
			Best_name = string("SA_BESTSOL_",g,".csv")
			Best_sol = open(Best_name, "w") #makes the CSV with above name for each step
			write(Best_sol, "Depot Selected, Solution, Temp, Step, GHG \r\n")
			for i in 2:length(flow_y)
				write(Best_sol,"$(flow_y[i]), $Best, $m, $g, $(e_step[g]) \r\n")
			end
			close(Best_sol)
		end
	
	#limit to number of unimproved ones
	if count == count_limit
		print("COUNT REACHED, BREAK!\r\n")
		break
	end

	end
end
end

total_time = time() - start
write(Solutions,"Total Time:, $(total_time) \r\n")
close(Solutions)
#close(SA_bases)
#close(swaps)
