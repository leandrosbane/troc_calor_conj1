use rfluids::prelude::*;
use std::error::Error;

const T_MIN: f64 = 250.0; // K
const T_MAX: f64 = 750.0; // K
const N: usize = 201;

#[derive(Clone, Copy, Debug)]
pub struct PropriedadesAr {
    pub rho: f64,
    pub cp: f64,
    pub k: f64,
    pub mu: f64,
}

pub struct TabelaAr {
    valores: Vec<PropriedadesAr>,
    passo: f64,
}

impl TabelaAr {
    // Executada uma única vez na inicialização
    pub fn nova(pressao: f64) -> Result<Self, Box<dyn Error>> {
        let passo = (T_MAX - T_MIN) / (N - 1) as f64;

        let mut ar = Fluid::from(Pure::Air).in_state(
            FluidInput::pressure(pressao),
            FluidInput::temperature(T_MIN),
        )?;

        let mut valores = Vec::with_capacity(N);

        for i in 0..N {
            let t = T_MIN + i as f64 * passo;

            ar.update(FluidInput::pressure(pressao), FluidInput::temperature(t))?;

            valores.push(PropriedadesAr {
                rho: ar.density()?,
                cp: ar.specific_heat()?,
                k: ar.conductivity()?,
                mu: ar.dynamic_viscosity()?,
            });
        }

        Ok(Self { valores, passo })
    }

    // Consulta rápida, sem chamar o CoolProp
    // Recebe temperatura em graus Celsius
    pub fn consultar(&self, temperatura: f64) -> PropriedadesAr {
        let t = temperatura + 273.15;

        assert!(
            t.is_finite() && (T_MIN..=T_MAX).contains(&t),
            "Temperatura fora dos limites da tabela: {temperatura} °C"
        );

        let indice = (((t - T_MIN) / self.passo) as usize).min(self.valores.len() - 2);

        let t0 = T_MIN + indice as f64 * self.passo;
        let f = (t - t0) / self.passo;

        let a = self.valores[indice];
        let b = self.valores[indice + 1];

        PropriedadesAr {
            rho: a.rho + f * (b.rho - a.rho),
            cp: a.cp + f * (b.cp - a.cp),
            k: a.k + f * (b.k - a.k),
            mu: a.mu + f * (b.mu - a.mu),
        }
    }
}
