plugins {
    kotlin("multiplatform") version "1.5.21"
}

group = "cc.eroute"
version = "1.0-SNAPSHOT"

repositories {
    mavenCentral()
}

kotlin {
    wasm32 {
        binaries {
            executable {
                entryPoint = "lemmingapex.trilateration.main"
            }
        }
    }
    jvm {
        compilations.all {
            kotlinOptions.jvmTarget = "1.8"
        }
        testRuns["test"].executionTask.configure {
            useJUnit()
        }
    }
    js(IR) {
        binaries.executable()
        browser()
    }
    val hostOs = System.getProperty("os.name")
    val isMingwX64 = hostOs.startsWith("Windows")
    val nativeTarget = when {
        hostOs == "Mac OS X" -> macosX64("native"){
            binaries {
                executable{
                    entryPoint = "trilateration.main"
                }
            }
        }
        hostOs == "Linux" -> linuxX64("native"){
            binaries {
                executable{
                    entryPoint = "trilateration.main"
                }
            }
        }
        isMingwX64 -> mingwX64("native"){
            binaries {
                executable{
                    entryPoint = "trilateration.main"
                }
            }
        }
        else -> throw GradleException("Host OS is not supported in Kotlin/Native.")
    }

    
    sourceSets {
        val commonMain by getting
        val commonTest by getting {
            dependencies {
                implementation(kotlin("test"))
            }
        }
        val jvmMain by getting
        val jvmTest by getting
        val jsMain by getting
        val jsTest by getting
        val nativeMain by getting
        val nativeTest by getting
        val wasm32Main by getting{
            dependencies {
                implementation(kotlin("stdlib"))
                implementation("org.jetbrains.kotlinx:kotlinx-serialization-runtime-wasm32:1.0-M1-1.4.0-rc")
            }
        }

    }
}
