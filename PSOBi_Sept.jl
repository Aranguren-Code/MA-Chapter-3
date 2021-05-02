# Loading the necessary packages
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
using CSV
using JuMP
using CPLEX
using DataFrames
using DelimitedFiles
#-------------------------------------------------------------------------------
# IF FIRST RUN, CHANGE TO 1
first = 1
if first == 1 
	print("Loading Model \r\n")
	include("getdata.jl");
	include("Load_Bi.jl");
else 
	include("reset_SA.jl");
end
#-------------------------------------------------------------------------------

Solutions = open("PSO Solutions.csv", "w")
PSO_moves = open("PSO_moves.csv", "w")
write(Solutions,"step, ite, Depot, YFLOW , next_energy, delta_e, Best, GHG \r\n")
start = time() 

write(PSO_moves, "ite,")
for i in 1:D
	write(PSO_moves, "$i,")
end
write(PSO_moves, "\r\n")#next row

current_energy = (maximum(ShortCost) * sum(Demand)) #third party supplier cost 
Best = current_energy #initial solution
next_energy = 0 #initialize variable
flag = 1
count = 0
#-----------------------------------------------------------
#-----------------------------------------------------------
# PSO Parameters 
#Number of Partibles
#-----------------------------------------------------------
# parts = 10
factor = 0.1 # for later make function?
parts = Int(ceil(D*factor))
c1 = 0.9 #local
c2 = 1-c1 #global
w = 1 #inertia
# If Best not found within n iterations, then stop :) 
n = 10
#------------------------------------------------------------
#-----------------------------------------------------------

mutable struct particles
	loc::Int
	locbest::Int
	supply::Float64
    lbest::Float64
end

base = Array{particles,1}(undef, parts)

#Starting all particles at 0
for i in 1:parts
    x = particles(0,0,0.0,current_energy)
    base[i] = x
end

#-----------------------------------------------------------------------
m=1 
	
for g in 1:20 #emissions steps
global count = 0
#Change Emission Limit----------------------------------------------
JuMP.fix(Epoint, e_step[g])
print("Emission step $g for $(e_step[g])  \r\n")
#-------------------------------------------------------------------

while count <= n
	if base[1].loc == 0 #if non selected, initialize randomly
		print("Initial selection: \r\n")
		for i in 1:parts
			r = ceil(D*rand())
			global x = particles(r,0,0.0,current_energy)
			global base[i] = x
			print(base[i].loc)
			print(", ")
			ID = base[i].loc
			JuMP.fix(D_Model[ID], 1) #Change model
		end
	else #velocity and change of 
		for i in 1:parts
			r1 = rand()
			r2 = rand()
				if base[i].supply>0
					print("Supply above zero!")
					vp = Int(ceil((c1*r1*(base[i].locbest - base[i].loc)))) #how to add global	
					if vp > 0
						vp = w + vp
					else #if the vp is negative
						vp = - w - vp
					end

					locp = base[i].loc + vp #new location up
					print("Changed $(base[i].loc), with $vp, to $locp")
				else
					j = i
					found = 0
					while found == 0
						if j == parts 
							break
						end
						if base[j].locbest > 0 
							found = 1
						end
						j+=1
					end
					print("Supply is zero, checked sorroundings, found $j.\r\n")
					vp = Int(w + ceil(c1*r1*(base[i].locbest - base[j].loc)))

					if vp > 0
						vp = w + vp
					else
						vp = - w - vp
					end
					locp = base[j].loc + vp
					print("Change $j, with $vp, to $locp\r\n")
				end			
				
			#Check if new location is out of bounds
			if locp > D #if outside of bounds
				print("Outside of bounds, $locp.")
				diff = locp - D
				locp = diff #we want to go back to beginning instead
				print("Fixed, new location, $locp. \r\n")
			elseif locp < 0
				locp = D + locp #go back from end, since locp is negative
			elseif locp == 0 # if locp is 0 
				locp = Int(ceil(D*rand())) #Select a random location
			end
			global base[i].loc = locp #update location
			JuMP.fix(D_Model[locp], 1)#Change model
		end
	end

	#-------------------------------------------------------------------------------------

	optimize!(HS_P3)

	#-----------------------------------------------------------------------------------
	# Energies and acceptance criteria

	global next_energy = JuMP.objective_value(HS_P3)
	global delta_e = current_energy - next_energy
	print("Delta_e: ")
	print(delta_e)
	print("\r\n")
	
	#TO Keep track of moves, may remove later
	#-----------------------------------------------------------------------------------#
	#to record selection
	b = zeros(Int,parts)
	for p in 1:parts
		global b[p] = base[p].loc
	end
		
	write(PSO_moves, "$m,")

	for i in 1:D
		if i in b
			write(PSO_moves, "1,")
		else
			write(PSO_moves, "0,")	
		end
	end

	write(PSO_moves, "$next_energy \r\n")
	
# #-----------------------------------------------------------------------------------------#
# #							 Accepting criteria 
# #-----------------------------------------------------------------------------------------#
	if delta_e >= 0
		global flag = 1
		global current_energy = next_energy
#--------------------------------------------------------
		#base = z_base #zero out base
		y_avg = 0
		#flow_y = zeros(Int64,1)
		#write(PSO_bases, "ite,particle, location, locbest, supply, global_best \r\n" )
		for p in 1:parts
			ID = base[p].loc
				for k in 1:C
					for s in 1:S
						if JuMP.value(FLOWY[ID,k,s]) > 0
							y_avg += JuMP.value(FLOWY[ID,k,s])*prob[s]
						end
					end
				end
				print(" Particle $p, Depots $(base[p].loc) ($(JuMP.value(D_Model[base[p].loc]))) has a supply of $y_avg \r\n")
				if y_avg > 0 
						print(ID)
						print(", ")
						print(y_avg)
						print("\r\n")
						
					if base[p].supply <= y_avg 
						print("Update best location: Update $(base[p].locbest) to $(ID) \r\n")
						global base[p].supply = y_avg #update best supply 
						global base[p].locbest = ID
						global base[p].lbest = current_energy
						#write(PSO_bases, "$m, $p, $(base[p].loc),$(base[p].locbest), $(base[p].supply), $(Best)\r\n")
					elseif current_energy < Best #if there improvement in investment
						global base[p].locbest = ID
						global base[p].lbest = current_energy #keep track of investment
					end
				else 
					print("Particle $p reset \r\n")
					global base[p].supply = 0
				end
			y_avg = 0
		end

		for p in 1:parts
			if base[p].supply == 0 
				ID = base[p].loc #move away more because of zero supply, use base[i].supply???
				print("Model reset depots: particle $p, Depot $(base[p].loc) \r\n")
				JuMP.fix(D_Model[ID],0)
			end
		end	

	elseif delta_e < 0 
		global flag = 0 # add to velocity??? 
		#Nothing is updated
		print("NO IMPROVEMENT \r\n")
		include("reset_SA.jl")
	end

	#write(Solutions,"ite, next_energy, delta_e, Best\r\n")
	if next_energy < Best 
		global Best = next_energy #update best
		for p in 1:parts
			if base[p].supply > 0 
				Depot_ID = base[p].loc 
				write(Solutions, "$g, $m, $(Depot_ID), $(base[p].supply),$next_energy, $delta_e, $Best, $(e_step[g])\r\n")
			end
		end
		count = 0
	end

	global m+=1 #keep track of iterations
	global count += 1

end

end

 total_time = time() - start
 write(PSO_moves,"Total Time:, $(total_time) \r\n")
 close(Solutions)
 close(PSO_moves)
