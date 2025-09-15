use std::fs;
use std::collections::HashMap;

pub struct ParameterReader
{
    content: String,
    // parameters: Vec<String>,
    parameters: &'static[&'static str]
}

// fn build_new_map(names: &[&'static str], content: &String) -> Result<HashMap<&'static str, String>,String>
// {
//     let mut parameter_map: HashMap<&'static str, String> = HashMap::new();
//     for line in content.lines()
//     {
//         for name in names
//         {
//             if line.contains(name)
//             {
//                 let Some((_, value)) = line.split_once(":") else
//                 {
//                     return Err(String::from("Bad delimiter"));
//                 };
//                 let value = value.trim();
//                 parameter_map.insert(name, value.to_owned());
//             }
//         }
//     }
//     for name in names
//     {
//         if !parameter_map.contains_key(name)
//         {
//             return Err(format!("Missing parameter: {name}"));
//         };
//     }
//     Ok(parameter_map)
// }

// pub fn read_parameters(args: &[String], parameters: &[&'static str]) -> Result<HashMap<&'static str, String>, String>
// {
//     if args.len() < 2
//     {
//         return Err(String::from("Not enoguh arguments: need filename containing parameters"));
//     }
//     let file_name = &args[1];
//     let content = match fs::read_to_string(file_name)
//     {
//         Ok(value) => value,
//         Err(err)  => return Err(err.to_string())
//     };

//     let parameter_dict = build_new_map(parameters, &content)?;

//     Ok(parameter_dict)
// }

impl ParameterReader
{
    pub fn build(args: &[String], parameters: &'static[&'static str]) -> Result<Self, String>
    {
        if args.len() < 2
        {
            return Err(format!("Not enough arguments: Usage: {} parameter_file", args[0]));
        }
        let file_name = &args[1];
        let content   = match fs::read_to_string(file_name)
        {
            Ok(value) => value,
            Err(err)  => return Err(err.to_string())
        };
        // let parameters: Vec<String> = parameters.iter().map(|str| str.to_string()).collect();
        
        let reader = Self {content, parameters};
        Ok(reader)
    }

    pub fn parse_parameters(&self, delimiter: &'static str) -> Result<HashMap<&'static str, String>, String>
    {
        let mut parameter_map: HashMap<&'static str, String> = HashMap::new();

        for line in self.content.lines()
        {
            for name in self.parameters
            {
                if line.contains(name)
                {
                    let Some((_, value)) = line.split_once(delimiter) else
                    {
                        return Err(String::from("Bad delimiter"));
                    };
                    let value = value.trim();
                    parameter_map.insert(name, value.to_owned());
                }
            }
        }

        for name in self.parameters
        {
            if !parameter_map.contains_key(name)
            {
                return Err(format!("Missing parameter: {name}"));
            };
        }
        Ok(parameter_map)
    }
}
 