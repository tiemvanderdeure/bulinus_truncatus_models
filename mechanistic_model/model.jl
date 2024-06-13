using LambertW

# Basic functions
quadratic(t, T0, Tm, q) = q * (t - T0) * (Tm - t)
gdd_function(t, T0, GDD) = ifelse(t > T0, (t-T0)/GDD, Float32(0))

# Define pop growth rate function
function get_pop_growth_rate!(output, ts = temperatures, traits_values = traits_values_f32)
    let 
        egglaying_T0 = traits_values["egglaying_T0"] 
        egglaying_Tm = traits_values["egglaying_Tm"]
        egglaying_q = traits_values["egglaying_q"]
        hatching_T0 = traits_values["hatchingtime_T0"]
        hatching_GDD = traits_values["hatchingtime_GDD"]
        maturation_T0 = traits_values["maturation_T0"]
        maturation_GDD = traits_values["maturation_GDD"]
        lifespan_T0 = traits_values["lifespan_T0"]
        lifespan_Tm = traits_values["lifespan_Tm"]
        lifespan_q = traits_values["lifespan_q"]
        hatching_success = traits_values["hatchingsuccess_q"]
        egg_death_rate = traits_values["egg_death_rate"]

        map!(output, ts) do t
            if ismissing(t) typeof(t)(NaN)
            else
                DR = (Float32(1) / exp(quadratic(t, lifespan_T0, lifespan_Tm, lifespan_q)))
                # if there is no egglaying, just return the death rate or hatching or maturing
                if (t < egglaying_T0) || (t > egglaying_Tm) || t < hatching_T0 || t < maturation_T0 
                    return - DR
                else 
                    ELR = quadratic(t, egglaying_T0, egglaying_Tm, egglaying_q)
                    hatching_time = Float32(1) / gdd_function(t, hatching_T0, hatching_GDD)/7
                    mat_time = Float32(1) / gdd_function(t, maturation_T0, maturation_GDD)
                    
                    # proportion of eggs that survives to maturity
                    surv_to_mature = hatching_success * exp.(-(mat_time * DR + hatching_time * egg_death_rate))
                    age_at_maturity = mat_time + hatching_time # avg time for an egg to reach maturity

                    if surv_to_mature > 0 && exp(age_at_maturity * DR) < Inf
                        return lambertw(ELR * age_at_maturity * surv_to_mature * exp(age_at_maturity * DR)) / age_at_maturity - DR
                    else
                        return - DR
                    end
                end
            end
        end
    end

    return
end 

function get_pop_growth_rate(ts = temperatures, traits_values = traits_values_f32)
    output = similar(ts)
    get_pop_growth_rate!(output, ts, traits_values)
    return output
end


