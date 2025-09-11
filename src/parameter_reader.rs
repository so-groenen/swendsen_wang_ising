use std::fs;
use std::collections::HashMap;
 

fn build_new_map(names: &[&'static str], content: &String) -> Result<HashMap<&'static str, String>,String>
{
    let mut parameter_map: HashMap<&'static str, String> = HashMap::new();
    for line in content.lines()
    {
        for name in names
        {
            if line.contains(name)
            {
                let Some((_, value)) = line.split_once(":") else
                {
                    return Err(String::from("Bad delimiter"));
                };
                let value = value.trim();
                parameter_map.insert(name, value.to_owned());
            }
        }
    }
    for name in names
    {
        if !parameter_map.contains_key(name)
        {
            return Err(format!("Missing parameter: {name}"));
        };
    }
    Ok(parameter_map)
}

pub fn read_parameters(args: &[String], parameters: &[&'static str]) -> Result<HashMap<&'static str, String>, String>
{
    if args.len() < 2
    {
        return Err(String::from("Not enoguh arguments: need filename containing parameters"));
    }
    let file_name = &args[1];
    let content = match fs::read_to_string(file_name)
    {
        Ok(value) => value,
        Err(err)  => return Err(err.to_string())
    };

    let parameter_dict = build_new_map(parameters, &content)?;

    Ok(parameter_dict)
}


 