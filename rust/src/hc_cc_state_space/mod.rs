// Copyright (c) 2017 - for information on the respective copyright
// owner see the NOTICE file and/or the repository
//
//     https://github.com/hbanzhaf/steering_functions.git
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

pub mod configuration;
pub mod hc_cc_circle;

pub use configuration::{
    configuration_aligned, configuration_distance, configuration_equal, Configuration,
};
pub use hc_cc_circle::{
    center_distance, configuration_on_hc_cc_circle, HcCcCircle, HcCcCircleParam,
};
