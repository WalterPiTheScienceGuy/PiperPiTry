// Arduino for fiber photometry
// To be controlled by python program

const int SHOCKER = 11; // shock US on pin 11
const int AUDIO = 9; // audio CS on pin 9
const int FIBER = 5; // fiber photometry event signal on pin 5
const int IR = 3; // infrared LED cue on pin 3
int incomingByte;      // a variable to read incoming serial data into

void setup() {
  // put your setup code here, to run once:
  Serial.begin(9600);
  pinMode(SHOCKER, OUTPUT); // set pin to digital output
  pinMode(AUDIO, OUTPUT); // set pin to digital output
  pinMode(FIBER, OUTPUT); 
  pinMode(IR, OUTPUT);
}

void loop() {
  // put your main code here, to run repeatedly:
  // see if there's incoming serial data:
  if (Serial.available() > 0) {
    // read the oldest byte in the serial buffer:
    incomingByte = Serial.read();
    // if it's a capital A (ASCII 65), turn on the LED:
    if (incomingByte == 'A') {
      digitalWrite(IR, HIGH);
    }
    // if it's a B (ASCII 66) turn off the LED:
    if (incomingByte == 'B') {
      digitalWrite(IR, LOW);
    }
    // if it's a capital C (ASCII 67), turn on the fiber pulse:
    if (incomingByte == 'C') {
      digitalWrite(FIBER, HIGH);
    }
    // if it's a D (ASCII 68) turn off the fiber pulse:
    if (incomingByte == 'D') {
      digitalWrite(FIBER, LOW);
    }
    // if it's a capital E (ASCII 69), turn on the AUDIO:
    if (incomingByte == 'E') {
      digitalWrite(AUDIO, HIGH);
    }
    // if it's a F (ASCII 70) turn off the AUDIO:
    if (incomingByte == 'F') {
      digitalWrite(AUDIO, LOW);
    }
    // if it's a capital G (ASCII 71), turn on the SHOCKER:
    if (incomingByte == 'G') {
      digitalWrite(SHOCKER, HIGH);
    }
    // if it's a H (ASCII 72) turn off the SHOCKER:
    if (incomingByte == 'H') {
      digitalWrite(SHOCKER, LOW);
    }
  }
}
