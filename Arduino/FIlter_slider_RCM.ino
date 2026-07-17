#include <SoftwareSerial.h>

int UV_fwd_pin = 34;
int UV_bwd_pin = 35;
int UV_state = 36;
int UV_RCM = 50;

int blue_fwd_pin = 38;
int blue_bwd_pin = 39;
int blue_state = 37;
int blue_RCM = 51;

int orange_fwd_pin = 40;
int orange_bwd_pin = 41;
int orange_state = 42;
int orange_RCM = 52;

int quad_fwd_pin = 44;
int quad_bwd_pin = 45;
int quad_state = 43;

int UV_signal = 0;
int blue_signal = 0;
int orange_signal = 0;

void setup() {
  // put your setup code here, to run once:
  Serial.begin(9600);  // open connection to PC for debugging
  Serial.write("Hello");
  pinMode(UV_fwd_pin, OUTPUT);
  pinMode(UV_bwd_pin, OUTPUT);
  pinMode(UV_state, INPUT);
  pinMode(UV_RCM, INPUT);

  pinMode(blue_fwd_pin, OUTPUT);
  pinMode(blue_bwd_pin, OUTPUT);
  pinMode(blue_state, INPUT);
  pinMode(blue_RCM, INPUT);

  pinMode(orange_fwd_pin, OUTPUT);
  pinMode(orange_bwd_pin, OUTPUT);
  pinMode(orange_state, INPUT);
  pinMode(orange_RCM, INPUT);

  pinMode(quad_fwd_pin, OUTPUT);
  pinMode(quad_bwd_pin, OUTPUT);
  pinMode(quad_state, INPUT);

}

void loop() {
  // put your main code here, to run repeatedly:
  UV_signal = digitalRead(UV_RCM);
  blue_signal = digitalRead(blue_RCM);
  orange_signal = digitalRead(orange_RCM);
  //Serial.println(RCM_state);

  if (UV_signal == 1){
    digitalWrite(UV_bwd_pin, 1);
    digitalWrite(UV_fwd_pin, 0);

    digitalWrite(blue_bwd_pin, 0);
    digitalWrite(blue_fwd_pin, 1);
    
    digitalWrite(orange_bwd_pin, 0);
    digitalWrite(orange_fwd_pin, 1);

    digitalWrite(quad_bwd_pin, 0);
    digitalWrite(quad_fwd_pin, 1);

    Serial.println("ON");
  }

  else if (blue_signal == 1){
    digitalWrite(UV_bwd_pin, 0);
    digitalWrite(UV_fwd_pin, 1);

    digitalWrite(blue_bwd_pin, 1);
    digitalWrite(blue_fwd_pin, 0);
    
    digitalWrite(orange_bwd_pin, 0);
    digitalWrite(orange_fwd_pin, 1);

    digitalWrite(quad_bwd_pin, 0);
    digitalWrite(quad_fwd_pin, 1);

    Serial.println("ON");
  }

  else if (orange_signal == 1){
    digitalWrite(UV_bwd_pin, 0);
    digitalWrite(UV_fwd_pin, 1);

    digitalWrite(blue_bwd_pin, 0);
    digitalWrite(blue_fwd_pin, 1);
    
    digitalWrite(orange_bwd_pin, 1);
    digitalWrite(orange_fwd_pin, 0);

    digitalWrite(quad_bwd_pin, 0);
    digitalWrite(quad_fwd_pin, 1);

    Serial.println("ON");
  }
  else {
    digitalWrite(UV_bwd_pin, 0);
    digitalWrite(UV_fwd_pin, 1);

    digitalWrite(blue_bwd_pin, 0);
    digitalWrite(blue_fwd_pin, 1);
    
    digitalWrite(orange_bwd_pin, 0);
    digitalWrite(orange_fwd_pin, 1);

    digitalWrite(quad_bwd_pin, 1);
    digitalWrite(quad_fwd_pin, 0);
    Serial.println("quad");
  }
  
  delay(200); 
  



}
