using JuMP

for i in 1:D
    JuMP.fix(D_Model[i],0)
end
print(" \r\n")
print("Successfull! \r\n")
print("All Depots have been set to zero (variable D_Model[j]=0)")
print(" \r\n")